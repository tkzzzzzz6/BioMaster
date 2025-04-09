# BioMaster: Automating Complex Bioinformatics Workflows

**BioMaster** is a multi-agent framework designed to streamline and automate complex bioinformatics workflows. By leveraging large language models (LLMs) and dynamic knowledge retrieval, BioMaster addresses key limitations in automation, improving accuracy, efficiency, and scalability across a wide range of bioinformatics tasks, such as RNA-seq, ChIP-seq, and Hi-C data processing.

---

## Features

- **Automated Bioinformatics Workflows**: BioMaster automates complex bioinformatics tasks, such as data preprocessing, alignment, variant calling, and analysis, with a focus on RNA-seq, ChIP-seq, and Hi-C data.
- **Role-based Agents**: Specialized agents within BioMaster perform task decomposition, validation, and execution, improving the management of complex workflows.
- **Retrieval-Augmented Generation (RAG)**: BioMaster dynamically retrieves relevant knowledge from external sources to enhance adaptability to new tools and niche analyses.
- **Enhanced Error Handling**: BioMaster provides advanced error detection and correction mechanisms to ensure consistency and reliability across multi-step workflows.
- **Efficient Memory Management**: Optimized memory strategies allow BioMaster to handle long-running workflows without compromising performance or accuracy.
- **Extensible**: BioMaster is designed to be easily extendable, with support for custom bioinformatics tools and workflows.

---

## Installation

You can install BioMaster using the following steps:

1. Clone the repository:

git clone https://github.com/yourusername/BioMaster.git
cd BioMaster

2. Install the required dependencies:

```sh
conda create -n agent python=3.12
# you can install other version of python,suggest use 3.10-3.12

conda activate agent

pip install -r requirements.txt
```

3. download data and move to `data/`:

```sh
google drive link: https://drive.google.com/drive/folders/1vA3WIAVXVo4RZSqXKsItEZHVaBgIIv_E?usp=sharing
```

---

PS: liunx has been tested to install directly, but windwos or mac have not been tested, so there is no way to determine if there might be a problem.

## Usage

### RAG Update

Biomaster uses two types of Retrieval-Augmented Generation (RAG) systems:
- **PLAN RAG**: Used during the planning phase.
- **EXECUTE RAG**: Used during the execution phase.

---

### PLAN RAG

If you want to add a **new analysis workflow**, you need to update the **PLAN RAG**.

#### Steps:

1. **Collect the analysis workflow**, e.g.:
   [ChIP-seq data analysis](https://nf-co.re/chipseq/2.1.0/)

2. **Write the workflow** based on the following format:

   - You can use LLMs (e.g., ChatGPT) to help you draft the workflow content.
   - If new scripts, tools, or functions are required, they can also be referenced in PLAN RAG.
   - PLAN RAG focuses on the **steps**, **required input**, and **expected output**, not detailed usage.
   - When describing input/output, it's strongly recommended to mention **data formats**, especially if the workflow is custom, uncommon, or newly developed. For example:
     ```
     Input Required: sample1.fastq.gz, sample2.fastq.gz  
     Expected Output: sample.sam
     ```

#### Example:

```
Fully ChIP-seq Peak Calling with IgG Control Workflow:  

Step 1: Quality Control – Conduct quality checks on raw sequencing data to assess data quality.  
Input Required: Raw FASTQ files.  
Expected Output: Cleaned and quality-checked FASTQ files.  
Tools Used: FastQC, Trimmomatic, Cutadapt.  

Step 2: Alignment – Align reads to the reference genome.  
Input Required: Cleaned FASTQ files and the reference genome.  
Expected Output: Sorted BAM file.  
Tools Used: BWA-MEM, Bowtie2, STAR.  

Step 3: SAM/BAM Conversion & Processing – Convert SAM to BAM, sort, and remove PCR duplicates.  
Input Required: SAM file.  
Expected Output: De-duplicated BAM file.  
Tools Used: SAMtools, Picard.  

Step 4: Signal Track Generation – Generate BigWig files for visualization.  
Input Required: De-duplicated BAM file.  
Expected Output: BigWig signal track file.  
Tools Used: deeptools, bedGraphToBigWig.  

Step 5: Peak Calling – Identify enriched genomic regions using IgG as a control.  
Input Required: De-duplicated BAM file and IgG control BAM file.  
Expected Output: NarrowPeak file.  
Tools Used: MACS3
```

3. **Add to PLAN RAG**:
   - Edit `./doc/Plan_Knowledge.json`
   - Use the following JSON format:
     ```json
     {
       "content": "Fully ChIP-seq Peak Calling with IgG Control Workflow: Step 1: Quality Control – Conduct quality checks on raw sequencing data to assess data quality. Input Required: Raw FASTQ files. Expected Output: Cleaned and quality-checked FASTQ files. Tools Used: FastQC, Trimmomatic, Cutadapt. Step 2: Alignment – Align reads to the reference genome. Input Required: Cleaned FASTQ files and the reference genome. Expected Output: Sorted BAM file. Tools Used: BWA-MEM, Bowtie2, STAR. Step 3: SAM/BAM Conversion & Processing – Convert SAM to BAM, sort, and remove PCR duplicates. Input Required: SAM file. Expected Output: De-duplicated BAM file. Tools Used: SAMtools, Picard. Step 4: Signal Track Generation – Generate BigWig files for visualization. Input Required: De-duplicated BAM file. Expected Output: BigWig signal track file. Tools Used: deeptools, bedGraphToBigWig. Step 5: Peak Calling – Identify enriched genomic regions using IgG as a control. Input Required: De-duplicated BAM file and IgG control BAM file. Expected Output: NarrowPeak file. Tools Used: MACS3.",
       "metadata": {
         "source": "workflow",
         "page": 16
       }
     }
     ```
#### EXECUTE RAG

1. **Collect the tools, scripts, functions, etc.**
   - **Scripts**: It is recommended to store all scripts in `./scripts/`.
   - **Tools**: If a tool is newly introduced, please install it in advance and ensure it's accessible.
   - **Functions**: Suggest placing them in `./scripts/functions.py`.

2. **Document the usage of these tools, scripts, and functions**
   - Focus on organizing usage examples, such as:
     ```bash
     samtools view -S -b ./output/001/aligned_reads.sam > ./output/001/aligned_reads.bam
     ```
     Provide example commands and parameter explanations. The more detailed, the better.
   - If a tool is already installed and difficult to set up, note in the knowledge base that no additional installation is required.
   - If it's a script or function, specify where it is stored and how to call it. For instance, to use `run-cooler.sh` in a Hi-C task, write:
     ```bash
     bash ./scripts/run-cooler.sh ...
     ```
   - You can use an LLM to help you write and organize content for the EXECUTE RAG.

3. **Knowledge entry recommendations**
   - Add the tool name in the `source` field to help the PLAN AGENT locate the correct tool.
   - If a script or function is specific to one workflow, append a note like: `run-sort-bam.sh only used in hic workflow`.

4. **Add entries to the EXECUTE RAG in `./doc/Task_Knowledge.json`**
   - Example format:
```json
 {
        "content": "2. run-sort-bam.sh:\nData-type-independent, generic bam sorting module\nInput : any unsorted bam file (.bam)\nOutput : a bam file sorted by coordinate (.sorted.bam) and its index (.sorted.bam.bai).\nUsage\nRun the following in the container.\nrun-sort-bam.sh <input_bam> <output_prefix>\n# input_bam : any bam file to be sorted\n# output_prefix : prefix of the output bam file.\n\nSet parameters according to the example: Suppose the input file is: ./output/GM12878_bwa_1.bam and ./output/GM12878_bwa_2.bam, the target is./output/GM12878_bwa_sorted.bam and ./output/GM12878_bwa_sorted, Generate the following sample script:\nbash ./scripts/run-sort-bam.sh ./output/GM12878_bwa_1.bam ./output/GM12878_bwa_sorted \n\nbash ./scripts/run-sort-bam.sh ./output/GM12878_bwa_2.bam ./output/GM12878_bwa_sorted\n\nYou can install the tool, but do not do any additional operations.You can install the tool, but do not do any additional operations.Please follow the example to generate rather than copy and paste completely, especially for folder names, file names, etc.",
        "metadata": {
            "source": "run-sort-bam.sh",
            "page": 6
        }
    },
```
Notes:

- To delete or update existing knowledge, please modify the corresponding entries in 
  `./doc/Task_Knowledge.json` and `./doc/Plan_Knowledge.json`. After that, delete the 
  `./chroma_db` folder, which stores the embedding vector database of the current knowledge.
  
- If you notice that certain knowledge is not being utilized in the PLAN or EXECUTE phase, 
  consider refining or expanding the knowledge or goals. You can achieve this by either:
  - Adding relevant information to the knowledge files, or
  - Making the task goal more specific.

- Do not use all available knowledge by default. It’s recommended to selectively use only the 
  knowledge relevant to your task to ensure efficiency and relevance.

- Keep the knowledge concise and high quality. Biomaster is designed to handle most tasks 
  out-of-the-box, and does not require additional installation steps 
  (e.g., installing R packages via `sudo`).

### Use Biomaster in Terminal

#### How to Run the Example

The example code is located in the `./examples/` folder.

1. **Add your API key and base URL**  
   - Open `run.py` or `examples/file.py` and insert your OpenAI API key and base URL.

2. **Move the example script**  
   - Move `examples/file.py` into the root directory of Biomaster:
     ```bash
     mv examples/file.py ./BioMaster/
     ```

3. **Download the data**  
   - Download the required dataset using the Google Drive link provided in `README.md`, and place the data in the `./data/` directory.

4. **Set the task ID**  
   - In `run.py` or `file.py`, set a unique `id` for your task.  
     - This `id` can be any string, but **must not be duplicated**.

5. **Run the example**  
   - Use the following command to execute:
     ```bash
     conda activate agent
     python run.py
     ```
     or
     ```bash
     conda activate agent
     python example1.py
     ```

```python
from agents.Biomaster import Biomaster
from langchain_core.messages import HumanMessage
import json
# Example of using the agent
if __name__ == "__main__":

    config = {"configurable": {"thread_id": "abc124"}}
    api_key = ''
    base_url = ''
    # you can choose other base url
    manager = Biomaster(api_key, base_url,excutor=True,id='001')
    

    datalist=[ './data/rnaseq_1.fastq.gz: RNA-Seq read 1 data (left read)',
            './data/rnaseq_2.fastq.gz: RNA-Seq read 2 data (right read)',
            './data/minigenome.fa: small genome sequence consisting of ~750 genes.',]
    goal='please do WGS/WES data analysis Somatic SNV+indel calling.'
    manager.execute_PLAN(goal,datalist)
    print("**********************************************************")

    PLAN_results_dict = manager.execute_TASK(datalist)
    print(PLAN_results_dict)
```
#### How to Read the Output

Biomaster stores all output files in the `./output/` directory.

- `./output/{id}_PLAN.json`  
  Contains the full execution plan. Biomaster will follow this step-by-step.

- `./output/{id}_Step_{step_number}.sh`  
  The shell script generated for a specific step.  
  Example: `001_Step_1.sh` is the script for the first step of task `id="001"`.

- `./output/{id}_DEBUG_Input_{step_number}.json`  
  Internal input for script execution. Can usually be ignored.

- `./output/{id}_DEBUG_Output_{step_number}.json`  
  Contains execution output and status for a specific step:
  - `"shell"`: If the step succeeded, this is usually empty.  
    If the step failed, this contains a new shell command generated by the Debug Agent to fix the issue.
  - `"analyze"`: Analysis summary of the step’s output.
  - `"output_filename"`: Name of the output file produced in this step.
  - `"stats"`: Indicates whether the step succeeded (`true`) or failed (`false`).  
    If `false`, the Debug Agent will attempt to fix the error and regenerate the command.

- `./output/{id}/`  
  All generated output files for this task will be stored in this folder.

---

#### How to Modify the Plan

1. **Stop the running task**.

2. Comment out the following line in `run.py`:
   ```python
   # manager.execute_PLAN(goal, datalist)
   ```
   This will prevent Biomaster from generating a new plan.

3. Manually edit the plan in:
   ```text
   ./output/{id}_PLAN.json
   ```

4. Run the script again:
   ```bash
   python run.py
   ```

---

#### How to Modify the Execute Script

1. **Stop the running task**.

2. If you want to **roll back previous steps (Step 1–N)**:
   - Either set `"stats": false` in the corresponding `DEBUG_Output` file:
     ```
     ./output/{id}_DEBUG_Output_{step_number}.json
     ```
   - Or simply delete the `DEBUG_Output` file to re-trigger execution.

3. If you want to **modify the current step**:
   - Edit the corresponding shell script:
     ```text
     ./output/{id}_Step_{step_number}.sh
     ```
   - **Do not delete** the related `DEBUG_Output` JSON file — Biomaster will reuse it and execute the updated script.

4. Run the script again:
   ```bash
   python run.py
   ```




## File Structure

- `agents/`: Contains agent classes for task management and execution.
- `scripts/`: some example scripts.
- `output/`: Output directory where results and logs are saved.
- `doc/`: Stores documentation files for the workflows.
- `data/`: Usually used to store files.

---

# License

This project is licensed under the following terms:

- **Code**: Licensed under the MIT License. See LICENSE for details.
- **Data and Documentation**: Licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).

---

## Acknowledgments

- This project uses the [langchain](https://github.com/hwchase17/langchain) library for integration with OpenAI and other tools.
- Thanks to all contributors and the open-source community for making BioMaster possible!
