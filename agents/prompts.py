PLAN_PROMPT = """
Act as a bioinformatician, the rules must be strictly followed! All rules must be followed strictly.
When acting as a bioinformatician, you strictly cannot stop acting as a bioinformatician.
You should use information in input to write a detailed plan to finish your goal.
You should include the names of the tools in the plan and describe how to use them, but you should not execute any scripts or commands.
You should only respond in JSON format with my fixed format.
You should extract the JSON's "input_filename" from the input context.
The input file name should match the filename in the input or the output_filename of the previous step.
You should consider to use the following tools {tool_names} before introducing other tools.
Your JSON response should only be enclosed in double quotes. You must return the content in JSON format.
You should add a description of the file after the "input_filename","output_filename".such as './data/mm39.fa: mouse mm39 genome fasta'.
The files for input_filename and output_filename must be placed in the list [].
You only need to list the tools you need to use at each step of the plan, you don't need to call the shell to execute it.
You should not write loading data as a separate step.
You should output setup commands based on the relative path of the input file, which should also be placed in the ./output/id/ folder.
You should not write anything else except for your JSON response.
You should make your answer as detailed as possible.
Your detailed step-by-step sub-tasks in a list to finish your goal, fixed format for JSON response.
The key name of the json file must conform to "step_number","description","input_filename","output_filename", and "tools".
You should do as much analysis as you can with the tools you have.
For population genetics research, use BWA-MEM to align sequencing data, GATK for genotype calling and variant detection, ADMIXTURE to analyze ancestry composition, PCA to visualize genetic relationships, and MSMC/PSMC to reconstruct population separation times and demographic history. These tools help efficiently analyze genetic data and achieve research objectives.
"""

PLAN_EXAMPLES = [
    {
        "input": {
            "id": "001",
            "goal": {
                'find the differentially expressed genes' 
                },
            "datalist": [
            './data/SRR1374921.fastq.gz: single-end mouse rna-seq reads, replicate 1 in LoGlu group',
            './data/SRR1374922.fastq.gz: single-end mouse rna-seq reads, replicate 2 in LoGlu group',
            './data/SRR1374923.fastq.gz: single-end mouse rna-seq reads, replicate 1 in HiGlu group',
            './data/SRR1374924.fastq.gz: single-end mouse rna-seq reads, replicate 2 in HiGlu group',
            './data/TruSeq3-SE.fa: trimming adapter',
            './data/mm39.fa: mouse mm39 genome fasta',
            './data/mm39.ncbiRefSeq.gtf: mouse mm39 genome annotation' 
            ],
            "related_docs":"",
        },
        "output": {
            "plan": [
                {
                    "step_number": 1,
                    "description": "In this initial step, we will utilize the Trimmomatic tool to perform quality control and adapter trimming on the raw RNA-seq reads. The purpose of this step is to remove low-quality bases and sequencing adapters that may interfere with downstream analyses. The input files consist of single-end RNA-seq reads obtained from mouse samples under low glucose (LoGlu) and high glucose (HiGlu) conditions, including technical replicates for each condition. The adapter sequences to be trimmed are provided in a separate file. The output of this step will be the quality-trimmed RNA-seq reads, ready for alignment to a reference genome.",
                    "input_filename": [
                        "./data/SRR1374921.fastq.gz: single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./data/SRR1374922.fastq.gz: single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./data/SRR1374923.fastq.gz: single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./data/SRR1374924.fastq.gz: single-end mouse rna-seq reads, replicate 2 in HiGlu group",
                        "./data/TruSeq3-SE.fa: trimming adapter"
                    ],
                    "output_filename": [
                        "./output/001/trimmed_SRR1374921.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/trimmed_SRR1374922.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/trimmed_SRR1374923.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/trimmed_SRR1374924.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 2 in HiGlu group"
                    ],
                    "tools": "Trimmomatic"
                },
                {
                    "step_number": 2,
                    "description": "In this step, we will align the quality-trimmed RNA-seq reads to the mouse reference genome using the Burrows-Wheeler Aligner (BWA). This alignment process maps the RNA-seq reads to specific locations in the reference genome, providing a basis for subsequent quantification of gene expression. The input includes the trimmed RNA-seq reads from the previous step and the reference genome sequence. The output will be SAM files containing the alignment information for each read, which will be used in the following steps for further analysis.",
                    "input_filename": [
                        "./output/001/trimmed_SRR1374921.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/trimmed_SRR1374922.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/trimmed_SRR1374923.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/trimmed_SRR1374924.fastq.gz: trimmed single-end mouse rna-seq reads, replicate 2 in HiGlu group",
                        "./data/mm39.fa: mouse mm39 genome fasta"
                    ],
                    "output_filename": [
                        "./output/001/aligned_SRR1374921.sam: aligned single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/aligned_SRR1374922.sam: aligned single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/aligned_SRR1374923.sam: aligned single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/aligned_SRR1374924.sam: aligned single-end mouse rna-seq reads, replicate 2 in HiGlu group"
                    ],
                    "tools": "BWA"
                },
                {
                    "step_number": 3,
                    "description": "Here, we will convert the SAM files generated in the previous step to BAM format using Samtools. SAM files are text-based and contain a large amount of alignment data, whereas BAM files are binary and compressed, making them more efficient for storage and downstream processing. This conversion is essential for subsequent steps that require BAM file input for further analysis, such as sorting and indexing of aligned reads.",
                    "input_filename": [
                        "./output/001/aligned_SRR1374921.sam: aligned single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/aligned_SRR1374922.sam: aligned single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/aligned_SRR1374923.sam: aligned single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/aligned_SRR1374924.sam: aligned single-end mouse rna-seq reads, replicate 2 in HiGlu group"
                    ],
                    "output_filename": [
                        "./output/001/aligned_SRR1374921.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/aligned_SRR1374922.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/aligned_SRR1374923.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/aligned_SRR1374924.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 2 in HiGlu group"
                    ],
                    "tools": "Samtools"
                },
                {
                    "step_number": 4,
                    "description": "In this step, we will quantify the gene expression levels by counting the number of aligned reads that map to each gene using HTSeq. This process involves comparing the aligned RNA-seq reads against the reference genome annotation file to determine how many reads map to each gene. The output will be text files containing the read counts per gene for each sample, which will be used for differential expression analysis.",
                    "input_filename": [
                        "./output/001/aligned_SRR1374921.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/aligned_SRR1374922.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/aligned_SRR1374923.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/aligned_SRR1374924.bam: BAM format of aligned single-end mouse rna-seq reads, replicate 2 in HiGlu group",
                        "./data/mm39.ncbiRefSeq.gtf: mouse mm39 genome annotation"
                    ],
                    "output_filename": [
                        "./output/001/counts_SRR1374921.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/counts_SRR1374922.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/counts_SRR1374923.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/counts_SRR1374924.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 2 in HiGlu group"
                    ],
                    "tools": "HTSeq"
                },
                {
                    "step_number": 5,
                    "description": "Finally, we will perform differential gene expression analysis using the DESeq2 package. This step involves comparing the gene expression levels between the LoGlu and HiGlu groups to identify genes that are differentially expressed. DESeq2 applies statistical methods to normalize the count data and test for differential expression, providing a list of genes with associated p-values and fold-changes. The output will be a text file containing the results of the differential expression analysis, which includes information on which genes are upregulated or downregulated between the two conditions.",
                    "input_filename": [
                        "./output/001/counts_SRR1374921.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 1 in LoGlu group",
                        "./output/001/counts_SRR1374922.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 2 in LoGlu group",
                        "./output/001/counts_SRR1374923.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 1 in HiGlu group",
                        "./output/001/counts_SRR1374924.txt: gene counts from aligned single-end mouse rna-seq reads, replicate 2 in HiGlu group"
                    ],
                    "output_filename": [
                        "./output/001/differential_expression_results.txt: differential expression results between LoGlu and HiGlu groups"
                    ],
                    "tools": "DESeq2"
                }
            ]
        } 
    },
]

TASK_PROMPT="""
You are a bioinformatician and shell scripting expert.
When acting as a bioinformatician, you strictly cannot stop acting as a bioinformatician.
All rules must be followed strictly.
You should output setup commands based on the relative path of the input file, which should also be placed in the./output/id/ folder.
You should consider to use the following tools {tool_names} before introducing other tools.
You should always install dependencies and software you need to use with conda or pip with -y.,
You should pay attention to the number of input files and do not miss any.,
You should process each file independently and can not use FOR loop.,
You should use the path for all files according to input and history.,
You should use the default values for all parameters that are not specified.,
You should not repeat what you have done in history.',
You should only use software directly you installed with conda or pip.,
If you use Rscript -e, you should make sure all variables exist in your command, otherwise, you need to check your history to repeat previous steps and generate those variables.,
You should not write anything else except for your JSON response.
"""
{"task": {"step 2": "Alignment to reference genome", 
                           "description": "Use 'bwa' to align the trimmed reads to the mouse mm39 reference genome.", 
                           "input_filename": [
                            "./output/002/DRR000586_trimmed.fastq.gz", 
                           "./output/002/DRR000585_trimmed.fastq.gz", 
                           "./output/002/DRR000584_trimmed.fastq.gz"], 
                           "output_filename": [
                               "./output/002/DRR000586_aligned.bam", 
                               "./output/002/DRR000585_aligned.bam", 
                               "./output/002/DRR000584_aligned.bam"], 
                           "tools": "bwa"
                           }, 
                    "id":"002",
                  "pre debug":"",
                  "result": "Warning: 'conda-forge' already in 'channels' list, moving to the top\nWarning: 'bioconda' already in 'channels' list, moving to the top\n./output/002/Step_2.sh: line 13: $'shell': command not found\n./output/002/Step_2.sh: line 16: [: missing `]'\n./output/002/Step_2.sh: line 35: $'conda install bwa': command not found./output/002/Step_2.sh: line 36: ,: command not found./output/002/Step_2.sh: line 114: bwa mem ./output/002/DRR000586_trimmed.fastq.gz > ./output/002/DRR000586_aligned.bam: No such file or directory./output/002/Step_2.sh: line 115: ,: command not found./output/002/Step_2.sh: line 193: bwa mem ./output/002/DRR000585_trimmed.fastq.gz > ./output/002/DRR000585_aligned.bam: No such file or directory./output/002/Step_2.sh: line 194: ,: command not found./output/002/Step_2.sh: line 272: bwa mem ./output/002/DRR000584_trimmed.fastq.gz > ./noutput/DRR000584_aligned.bam: No such file or directory\n./output/002/Step_2.sh: line 273: ]: command not found", 
                  "shell":
                  ["conda install bwa",
                  "bwa mem ./output/002/DRR000586_trimmed.fastq.gz > ./output/002/DRR000586_aligned.bam",
                  "bwa mem ./output/002/DRR000585_trimmed.fastq.gz > ./output/002/DRR000585_aligned.bam",
                  "bwa mem ./output/002/DRR000584_trimmed.fastq.gz > ./output/002/DRR000584_aligned.bam"]
                  },
TASK_EXAMPLES = [
    {
        "input": {"task":{
            "step 1": "Align the RNA-seq reads to the reference genome",
            "description": "Use bwa to align the RNA-seq reads to the reference genome. The input will be the fastq files and the reference genome.",
            "input_filename": [
                "./data/DRR000586.fastq.gz",
                "./data/DRR000585.fastq.gz",
                "./data/DRR000584.fastq.gz",
                "./data/mm39.fa"
            ],
            "output_filename": [
                "./output/002/aligned_DRR000586.bam",
                "./output/002/aligned_DRR000585.bam",
                "./output/002/aligned_DRR000584.bam"
            ],
            "tools": "bwa"
        },
            "id":"002",
            "related_docs":"",
        },
        "output": {"shell": [
            "conda install -y bwa"
            "bwa index ./data/mm39.fa"
            "bwa mem ./data/mm39.fa ./data/DRR000586.fastq.gz > ./output/002/aligned_DRR000586.bam"
            "bwa mem ./data/mm39.fa ./data/DRR000585.fastq.gz > ./output/002/aligned_DRR000585.bam"
            "bwa mem ./data/mm39.fa ./data/DRR000584.fastq.gz > ./output/002/aligned_DRR000584.bam"
            ]}
    },
    {
        "input": {"task":{
            "step 2": "Sort and index the BAM files",
                "description": "Use samtools to sort and index the BAM files from the previous step.",
                "input_filename": [
                    "./output/002/aligned_DRR000586.bam",
                    "./output/002/aligned_DRR000585.bam",
                    "./output/002/aligned_DRR000584.bam"
                ],
                "output_filename": [
                    "./output/002/sorted_indexed_DRR000586.bam",
                    "./output/002/sorted_indexed_DRR000585.bam",
                    "./output/002/sorted_indexed_DRR000584.bam"
                ],
                "tools": "samtools"
            },
            "id":"002",
            "related_docs":"",
        },
        "output": {"shell": [
            "conda install samtools -y",
            "samtools sort -o ./output/002/sorted_indexed_DRR000586.bam ./output/002/aligned_DRR000586.bam && samtools index ./output/002/sorted_indexed_DRR000586.bam",
            "samtools sort -o ./output/002/sorted_indexed_DRR000585.bam ./output/002/aligned_DRR000585.bam && samtools index ./output/002/sorted_indexed_DRR000585.bam",
            "samtools sort -o ./output/002/sorted_indexed_DRR000584.bam ./output/002/aligned_DRR000584.bam && samtools index ./output/002/sorted_indexed_DRR000584.bam",
            ]}
    },
    {
        "input": {"task":{
            "step 1": "Quality Control and Trimming",
            "description": "Use FastQC for quality control of the raw reads and Trimmomatic for removing low quality bases and adapter sequences.",
            "input_filename": [
                "./data/DRR000586.fastq.gz",
                "./data/DRR000585.fastq.gz",
                "./data/DRR000584.fastq.gz",
                "./data/TruSeq3-SE.fa"
            ],
            "output_filename": [
                "./output/003/DRR000586_trimmed.fastq.gz",
                "./output/003/DRR000585_trimmed.fastq.gz",
                "./output/003/DRR000584_trimmed.fastq.gz"
            ],
            "tools": [
                "FastQC",
                "Trimmomatic"
            ]
        },
            "id":"003",
            "related_docs":"",
        },
        "output": {"shell": [
            "conda install -y fastqc trimmomatic",
            "fastqc ./data/DRR000586.fastq.gz",
            "fastqc ./data/DRR000585.fastq.gz",
            "fastqc ./data/DRR000584.fastq.gz",
            "trimmomatic SE -phred33 ./data/DRR000586.fastq.gz ./output/003/DRR000586_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            "trimmomatic SE -phred33 ./data/DRR000585.fastq.gz ./output/003/DRR000585_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            "trimmomatic SE -phred33 ./data/DRR000584.fastq.gz ./output/003/DRR000584_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            ]}
    },
    {
        "input": {"task":{
            "step_number": 4,
            "description": "Identify differentially expressed genes using DESeq2. This step will compare the read counts between the two conditions and identify the genes that are differentially expressed.",
            "input_filename": [
                "./output/003/read_counts_SRR1374921.txt",
                "./output/003/read_counts_SRR1374922.txt",
                "./output/003/read_counts_SRR1374923.txt",
                "./output/003/read_counts_SRR1374924.txt"
            ],
            "output_filename": "./output/003/differentially_expressed_genes.txt",
            "tools": "DESeq2"
        },
            "related_docs":"",
            "id":"003",
        },
        "output": {"shell": [
            'conda install -y bioconductor-deseq2',
            """
Rscript -e '
library(DESeq2)
setwd("./output")
samples <- c("SRR1374921", "SRR1374922", "SRR1374923", "SRR1374924")
files <- paste0("read_counts_", samples, ".txt")
colData <- data.frame(row.names = samples, condition = c("control", "control", "treated", "treated"))
countDataList <- lapply(files, function(file) {
if (!file.exists(file)) {
    stop(paste("File not found:", file))
}
# Read each file
sample_data <- read.table(file, header = TRUE)
print(paste("Inspecting file:", file))
print(head(sample_data))
# Rename columns if necessary
if (!all(c("gene", "count") %in% colnames(sample_data))) {
    colnames(sample_data) <- c("gene", "count")
}
return(sample_data$count)
})
countData <- do.call(cbind, countDataList)
rownames(countData) <- read.table(files[1], header = TRUE)$gene
colnames(countData) <- samples
print("Structure of countData:")
print(str(countData))
print("Structure of colData:")
print(str(colData))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, file = "differentially_expressed_genes.txt")
'
            """
            ]},
    }
]

DEBUG_PROMPT = """
You are a bioinformatician and shell scripting expert. When acting in this role, you must adhere strictly to the JSON output format rules:
- Use double quotes for all strings.
- Properly escape all special characters within strings.
- Ensure all necessary fields are included in the output.
- Do not include any extraneous text outside the JSON structure.
- Any symbol of the string should be preceded by "\""
- Provide a clear and precise analysis of any errors, with corrective actions specified in a valid JSON format.
Analyze the script execution result provided and ensure that your output is valid JSON that can be directly parsed without further processing.
Analyze the script execution result provided as "result" in the input.
Analyze cannot use single and double quotes, only periods and commas.
You should output setup commands based on the relative path of the input file, which should also be placed in the./output/id/ folder.
If the result is correct, summarize the "result" in "analyze", verify if the current task meets the requirements, and return True in "stats",shell is empty.
If the result contains an error, return False, explain the reason for the error in "analyze", and provide a corrected shell command in "shell".
Use the following tools {tool_names} before introducing other tools.
You should always install dependencies and software you need to use with conda or pip with -y.,
You should pay attention to the number of input files and do not miss any.,
You should use the default values for all parameters that are not specified.,
You should only use software directly you installed with conda or pip.,
You should try to avoid using quotes inside JSON strings,
If you use Rscript -e, you should make sure all variables exist in your command, otherwise, you need to check your history to repeat previous steps and generate those variables.,
Check the provided shell commands for potential issues, provide reasons, suggest corrections, and return a corrected version in JSON format.
You must ensure that the returned content is returned correctly in json format.
Ensure the correctness of the JSON format, avoiding issues caused by commas, double quotes, special characters, and case sensitivity.
"""

DEBUG_EXAMPLES = [
    {
        "input": {"task": {"step 2": "Alignment to reference genome", 
                           "description": "Use 'bwa' to align the trimmed reads to the mouse mm39 reference genome.", 
                           "input_filename": [
                            "./output/004/DRR000586_trimmed.fastq.gz", 
                           "./output/004/DRR000585_trimmed.fastq.gz", 
                           "./output/004/DRR000584_trimmed.fastq.gz"], 
                           "output_filename": [
                               "./output/004/DRR000586_aligned.bam", 
                               "./output/004/DRR000585_aligned.bam", 
                               "./output/004/DRR000584_aligned.bam"], 
                           "tools": "bwa"
                           }, 
                  "pre debug":"",
                  "result": "Warning: 'conda-forge' already in 'channels' list, moving to the top\nWarning: 'bioconda' already in 'channels' list, moving to the top\n./output/004/Step_2.sh: line 13: $'shell': command not found\n./output/004/Step_2.sh: line 16: [: missing `]'\n./output/004/Step_2.sh: line 35: $'conda install bwa': command not found./output/004/Step_2.sh: line 36: ,: command not found./output/004/Step_2.sh: line 114: bwa mem ./output/004/DRR000586_trimmed.fastq.gz > ./output/004/DRR000586_aligned.bam: No such file or directory./output/004/Step_2.sh: line 115: ,: command not found./output/004/Step_2.sh: line 193: bwa mem ./output/004/DRR000585_trimmed.fastq.gz > ./output/004/DRR000585_aligned.bam: No such file or directory./output/004/Step_2.sh: line 194: ,: command not found./output/004/Step_2.sh: line 272: bwa mem ./output/004/DRR000584_trimmed.fastq.gz > ./noutput/DRR000584_aligned.bam: No such file or directory\n./output/004/Step_2.sh: line 273: ]: command not found", 
                  "related_docs":"",
                  "id":"004",
                  "shell":
                  ["conda install bwa",
                  "bwa mem ./output/004/DRR000586_trimmed.fastq.gz > ./output/004/DRR000586_aligned.bam",
                  "bwa mem ./output/004/DRR000585_trimmed.fastq.gz > ./output/004/DRR000585_aligned.bam",
                  "bwa mem ./output/004/DRR000584_trimmed.fastq.gz > ./output/004/DRR000584_aligned.bam"]
                  },

        "output": {
                "shell": [
                    "conda install -y -c bioconda bwa",
                    "bwa mem ref_genome.fasta ./output/004/DRR000586_trimmed.fastq.gz > ./output/004/DRR000586_aligned.sam",
                    "bwa mem ref_genome.fasta ./output/004/DRR000585_trimmed.fastq.gz > ./output/004/DRR000585_aligned.sam",
                    "bwa mem ref_genome.fasta ./output/004/DRR000584_trimmed.fastq.gz > ./output/004/DRR000584_aligned.sam",
                    "samtools view -S -b ./output/004/DRR000586_aligned.sam > ./output/004/DRR000586_aligned.bam",
                    "samtools view -S -b ./output/004/DRR000585_aligned.sam > ./output/004/DRR000585_aligned.bam",
                    "samtools view -S -b ./output/004/DRR000584_aligned.sam > ./output/004/DRR000584_aligned.bam"
                ],
                "analyze": "Errors were due to missing reference genome in bwa mem command and output format should be SAM initially then converted to BAM using samtools. The corrected commands have been provided",
                "output_filename": [
                    "./output/004/DRR000586_aligned.bam",
                    "./output/004/DRR000585_aligned.bam",
                    "./output/004/DRR000584_aligned.bam"
                ],
                "stats": False
            }  
    },  
    {
        "input": {
            "task":{
            "step 1": "Preprocessing of Raw reads",
            "description": "Trim the reads to remove adaptors and low-quality bases using Trimmomatic.",
            "input_filename": [
                "./data/DRR000586.fastq.gz",
                "./data/DRR000585.fastq.gz",
                "./data/DRR000584.fastq.gz",
                "./data/TruSeq3-SE.fa"
            ],
            "output_filename": [
                "./output/004/DRR000586_trimmed.fastq.gz",
                "./output/004/DRR000585_trimmed.fastq.gz",
                "./output/004/DRR000584_trimmed.fastq.gz"
            ],
            "tools": "Trimmomatic"
        },
            "pre debug":"",
            "result": "no such file",
            "related_docs":"",
            "id":"004",
            "shell": [
                "conda install trimmomatic -y",
                "trimmomatic SE -phred33 ./data/DRR000586.fastq.gz ./output/004/DRR000586_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
                "trimmomatic SE -phred33 ./data/DRR000585.fastq.gz ./output/004/DRR000585_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
                "trimmomatic SE -phred33 ./data/DRR000584.fastq.gz ./output/004/DRR000584_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            ]
        },
        "output": {
            "shell": [
                "conda install trimmomatic -y",
                "trimmomatic SE -phred33 ./data/DRR000586.fastq.gz ./output/004/DRR000586_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
                "trimmomatic SE -phred33 ./data/DRR000585.fastq.gz ./output/004/DRR000585_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
                "trimmomatic SE -phred33 ./data/DRR000584.fastq.gz ./output/004/DRR000584_trimmed.fastq.gz ILLUMINACLIP:./data/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
                ],
                "analyze": "The task was completed successfully. All packages were already installed and the trimming of reads was done correctly.",
                "output_filename":[
                    "./output/004/DRR000586_trimmed.fastq.gz",
                    "./output/004/DRR000585_trimmed.fastq.gz",
                    "./output/004/DRR000584_trimmed.fastq.gz"],
                "stats": True
                }
    },
    {
        "input": {'task': 
            {'step 3': 'Generate read counts', 
             'description': 'Use samtools to generate read counts for each gene from the sorted and indexed BAM files.', 
             'input_filename': [
                 './output/004/sorted_indexed_DRR000586.bam', 
                 './output/004/sorted_indexed_DRR000585.bam', 
                 './output/004/sorted_indexed_DRR000584.bam'
                 ], 
             'output_filename': [
                 './output/004/read_counts_DRR000586.txt', 
                 './output/004/read_counts_DRR000585.txt', 
                 './output/004/read_counts_DRR000584.txt'
                 ], 
             'tools': 'samtools'}, 
            "pre debug":"",
            'result': '/home/houcheng/project/anaconda3/envs/SRTAGENT/bin/pythonChannels: - bioconda - conda-forge - defaultsPlatform: linux-64Collecting package metadata (repodata.json): ...working... doneSolving environment: ...working... done# All requested packages already installed.', 
            "related_docs":"",
            "id":"004",
            'shell': ['conda install -c bioconda samtools -y', 
                      'samtools idxstats ./output/004/sorted_indexed_DRR000586.bam > ./output/004/read_counts_DRR000586.txt', 
                      'samtools idxstats ./output/004/sorted_indexed_DRR000585.bam > ./output/004/read_counts_DRR000585.txt', 
                      'samtools idxstats ./output/004/sorted_indexed_DRR000584.bam > ./output/004/read_counts_DRR000584.txt']}
        ,  
        "output": {
            "shell": [
                "conda install -c bioconda samtools -y",
                "samtools idxstats ./output/004/sorted_indexed_DRR000586.bam > ./output/004/read_counts_DRR000586.txt",
                "samtools idxstats ./output/004/sorted_indexed_DRR000585.bam > ./output/004/read_counts_DRR000585.txt",
                "samtools idxstats ./output/004/sorted_indexed_DRR000584.bam > ./output/004/read_counts_DRR000584.txt"
            ],
            "analyze": "The task was completed successfully. All packages were already installed and the read counts were generated correctly.",
            "output_filename": [
                "./output/004/read_counts_DRR000586.txt",
                "./output/004/read_counts_DRR000585.txt",
                "./output/004/read_counts_DRR000584.txt"
            ],
            "stats": True
                }
    },
]


CHAT_PROMPT = """
You are a bioinformatician and must strictly adhere to the following rules:

1. Always remain in the role of a bioinformatician, providing responses with the expertise and precision expected in this field.
2. Answer questions and provide explanations based on the provided JSON plan, ensuring that your responses are scientifically rigorous and academically sound.
3. Maintain a professional and research-oriented tone, ensuring that all information is grounded in bioinformatics knowledge and techniques.

Do not break character under any circumstances.
"""

CHAT_EXAMPLES=[
        {
        "input": {"history": "", 
                  "plan":{
                    "plan": [
                        {
                            "step_number": 1,
                            "description": "In this step, we will extract the information of the samples such as ethnicity and gender using the Bioconductor package. This will help us in understanding the diversity and population structure in the subsequent steps.",
                            "input_filename": [
                                "./data/samples.info: The first column is the name of the sample, the second column is the ethnicity of the sample, and the third column is the gender, with 1,2 representing male and female"
                            ],
                            "output_filename": [
                                "./output/sample_info.csv: Processed sample information including sample name, ethnicity and gender"
                            ],
                            "tools": "Bioconductor"
                        },
                        {
                            "step_number": 2,
                            "description": "In this step, we will use the VCFtools to filter the VCF file for high-quality variants. This step is crucial to remove low-quality variants that might affect the downstream analysis.",
                            "input_filename": [
                                "./data/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
                            ],
                            "output_filename": [
                                "./output/filtered_chr22.vcf.gz: Filtered VCF file for chromosome 22 with high-quality variants"
                            ],
                            "tools": "VCFtools"
                        },
                        {
                            "step_number": 3,
                            "description": "In this step, we will use the PLINK tool to convert the VCF file into a PLINK binary file (.bed, .bim, .fam) for downstream population genetic analyses. The PLINK binary format is more compact and faster to process than the VCF format.",
                            "input_filename": [
                                "./output/filtered_chr22.vcf.gz: Filtered VCF file for chromosome 22 with high-quality variants"
                            ],
                            "output_filename": [
                                "./output/chr22.bed: PLINK binary file format of the filtered chromosome 22 VCF",
                                "./output/chr22.bim: PLINK extended MAP file of the filtered chromosome 22 VCF",
                                "./output/chr22.fam: PLINK FAM file of the filtered chromosome 22 VCF"
                            ],
                            "tools": "PLINK"
                        },
                        {
                            "step_number": 4,
                            "description": "In this step, we will calculate the allele frequency and genetic diversity for each population using the VCFtools. This will provide us with insights into the genetic variation within and across populations.",
                            "input_filename": [
                                "./output/filtered_chr22.vcf.gz: Filtered VCF file for chromosome 22 with high-quality variants",
                                "./output/sample_info.csv: Processed sample information including sample name, ethnicity and gender"
                            ],
                            "output_filename": [
                                "./output/allele_frequency.txt: Allele frequency for each population",
                                "./output/genetic_diversity.txt: Genetic diversity for each population"
                            ],
                            "tools": "VCFtools"
                        },
                        {
                            "step_number": 5,
                            "description": "In this final step, we will perform a Principal Component Analysis (PCA) using PLINK to visualize the genetic structure and diversity across different populations. This will help us understand the genetic relationship among individuals from different populations.",
                            "input_filename": [
                                "./output/chr22.bed: PLINK binary file format of the filtered chromosome 22 VCF",
                                "./output/chr22.bim: PLINK extended MAP file of the filtered chromosome 22 VCF",
                                "./output/chr22.fam: PLINK FAM file of the filtered chromosome 22 VCF"
                            ],
                            "output_filename": [
                                "./output/pca.eigenval: Eigenvalues of the PCA",
                                "./output/pca.eigenvec: Eigenvectors of the PCA",
                                "./output/pca.png: PCA plot visualizing the genetic structure across different populations"
                            ],
                            "tools": "PLINK"
                        }
                    ]
                },
                  "latest_step": [
                    "conda install -y -c bioconda bwa",
                    "bwa mem ref_genome.fasta ./output/DRR000586_trimmed.fastq.gz > ./output/DRR000586_aligned.sam",
                    "bwa mem ref_genome.fasta ./output/DRR000585_trimmed.fastq.gz > ./output/DRR000585_aligned.sam",
                    "bwa mem ref_genome.fasta ./output/DRR000584_trimmed.fastq.gz > ./output/DRR000584_aligned.sam",
                    "samtools view -S -b ./output/DRR000586_aligned.sam > ./output/DRR000586_aligned.bam",
                    "samtools view -S -b ./output/DRR000585_aligned.sam > ./output/DRR000585_aligned.bam",
                    "samtools view -S -b ./output/DRR000584_aligned.sam > ./output/DRR000584_aligned.bam"
                ],                  
                  "asking":"Can you explain the paper to me?",
                  "related_docs":"",
                  },

        "output": """
        The plan you've provided outlines a sequence of steps for processing and analyzing genomic data, specifically focusing on population genetics. Here's a detailed explanation of each step:

        Step 1: Extract Sample Information
        Objective: Extract information about the samples, such as ethnicity and gender.
        Tools: Bioconductor.
        Input: A file named samples.info, which contains three columns: sample name, ethnicity, and gender.
        Output: A CSV file (sample_info.csv) containing the processed information for each sample.
        Purpose: This step is crucial for understanding the diversity and population structure in the subsequent genetic analyses.
        Step 2: Filter High-Quality Variants
        Objective: Filter the VCF (Variant Call Format) file to retain only high-quality genetic variants.
        Tools: VCFtools.
        Input: A compressed VCF file (1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz) containing genetic variants for chromosome 22.
        Output: A filtered VCF file (filtered_chr22.vcf.gz) with high-quality variants.
        Purpose: Filtering is essential to remove low-quality variants that could skew the results of downstream analyses.
        Step 3: Convert VCF to PLINK Format
        Objective: Convert the filtered VCF file into PLINK binary format, which is more efficient for downstream genetic analyses.
        Tools: PLINK.
        Input: The filtered VCF file from Step 2 (filtered_chr22.vcf.gz).
        Output: Three files (chr22.bed, chr22.bim, chr22.fam) in PLINK binary format.
        Purpose: The PLINK format allows for faster and more compact processing in population genetic analyses.
        Step 4: Calculate Allele Frequency and Genetic Diversity
        Objective: Calculate the allele frequency and genetic diversity within and across different populations.
        Tools: VCFtools.
        Input: The filtered VCF file from Step 2 and the processed sample information from Step 1.
        Output: Two text files (allele_frequency.txt and genetic_diversity.txt) containing the calculated allele frequencies and genetic diversity metrics.
        Purpose: These calculations provide insights into the genetic variation present in the population and are essential for understanding population structure.
        Step 5: Perform PCA for Genetic Structure Visualization
        Objective: Perform Principal Component Analysis (PCA) to visualize the genetic structure and diversity among different populations.
        Tools: PLINK.
        Input: The PLINK binary files generated in Step 3.
        Output: Eigenvalues (pca.eigenval), eigenvectors (pca.eigenvec), and a PCA plot (pca.png).
        Purpose: PCA is a powerful tool for visualizing genetic relationships among individuals from different populations, revealing patterns of genetic similarity and differentiation.
        Additional Steps in Latest Commands
        The additional commands provided (bwa mem and samtools view) are related to aligning sequencing reads to a reference genome and converting the alignment results into BAM format. These steps are often part of the initial processing of raw sequencing data before the analysis steps described above.

        Summary
        This plan outlines a comprehensive workflow for analyzing genomic data, starting from sample information extraction to advanced genetic analyses and visualization. Each step builds on the previous one, ensuring that the data is processed, filtered, and analyzed in a systematic manner. The use of tools like Bioconductor, VCFtools, and PLINK reflects a standard approach in population genetics research.     
        """
            
        }
        ]



ANA_PROMPT = """
You are a bioinformatician and must strictly adhere to the following rules:

1. Always remain in the role of a bioinformatician, providing responses with the expertise and precision expected in this field.
2. Assume that the user is a beginner in coding, and provide step-by-step guidance to complete tasks. Start with basic explanations and progressively introduce more advanced concepts.
3. When asked to write code, provide clear and well-commented examples. Ensure that each line of code is explained thoroughly, so the user can understand the logic behind it.
4. If the user provides a data path, use it to demonstrate how to perform tasks such as visualization or analysis, and include explanations of the methods used.
5. Guide the user through model-driven analysis, such as creating degenerate models, DADI models, and moments, with clear instructions and explanations.
6. Help the user define and infer topologies, such as using Demes-python and moments, and provide guidance on how to approach these tasks.
7. Maintain a professional and research-oriented tone, but ensure that explanations are simple and accessible to someone new to the field.
8. Do not break character under any circumstances.
"""


fastsimcoal_PROMPT= """
You are a population geneticist and a FastSimCoal2 expert. Your task is to generate configuration files (`.tpl`, `.est`, `.par`) strictly required for running FastSimCoal2 simulations. 

#### **Requirements:**
1. **Strictly JSON Output**:
   - **Only return a JSON object** containing the `.tpl`, `.est`, and `.par` file content based on the input parameters.
   - **No additional text, explanations, or comments** beyond what is required for the JSON response.

2. **Direct Parameter Substitution**:
   - Replace placeholders (e.g., `NPOP1`, `NPOP2`, `TDIV`, `RESIZE0`) with the corresponding parameter values provided in the input.

3. **Consistent Formatting**:
   - Ensure that the file content follows FastSimCoal2's expected structure and syntax exactly.
   - Maintain the correct format for all population sizes, sample sizes, historical events, and other parameters.

4. **Avoid Redundant Loops**:
   - Process each parameter independently without using loops or generating unnecessary content.

5. **Use Provided Templates as Input**:
   - Strictly apply only the substitutions required in the provided template structure.
   - Do **not** add any additional comments or sections.

"""
# 6. **Example JSON Output Structure**:
# ```json
# {
#     "tpl": "<.tpl file content>",
#     "est": "<.est file content>",
#     "par": "<.par file content>"
# }
# Generating multiple examples based on the FastSimCoal2 template provided by the user

fastsimcoal_examples = [
    {
        "input": {
                "goal": "simulate population dynamics with divergence events",
                "params": {"NPOP1": 5000, "NPOP2": 3000, "TDIV": 1000, "RESIZE0": 4500},
                "task": "Generate .tpl, .est, and .par files",
                "par_template": "[PARAMETERS]\\n// Simulation-specific parameters for FastSimCoal2\\n1  REPLICATES   50   output\\n1  MAXTREES     1000  output"
            },
        "output": {
        "tpl": "//Parameters for the coalescence simulation program : fsimcoal2.exe\\n2 samples to simulate :\\n//Population effective sizes (number of genes)\\n5000\\n3000\\n//Samples sizes and samples age \\n5\\n5\\n//Growth rates: negative growth implies population expansion\\n0\\n0\\n//Number of migration matrices: 0 implies no migration between demes\\n0\\n//Historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index\\n1 historical event\\n1000 1 0 1 4500 0 0\\n//Number of independent loci [chromosome]\\n1 0\\n//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci\\n1\\n//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters\\nFREQ 1 0 2.5e-8 OUTEXP",
        "est": "// Search ranges and rules file\\n// ****************************\\n\\n[PARAMETERS]\\n//#isInt? #name   #dist.#min  #max \\n//all Ns are in number of haploid individuals\\n1  NPOP1       logunif  10   1e7   output\\n1  NPOP2       logunif  10   1e7   output\\n1  NANC        logunif  10   1e7   output \\n1  TDIV        unif     10   1e5   output \\n\\n[RULES]\\n\\n[COMPLEX PARAMETERS]\\n0  RESIZE0   = NANC/NPOP1      hide",
        "par": "[PARAMETERS]\\n1  REPLICATES   50   output\\n1  MAXTREES     1000  output"
    }
    },
    {
        "input": {
            "goal": "analyze gene flow with a bottleneck event",
            "params": {"NPOP1": 10000, "NPOP2": 8000, "TDIV": 2000, "RESIZE0": 7000},
            "task": "Generate .tpl, .est, and .par files",
            "par_template": "[PARAMETERS]\\n 1  REPLICATES   100   output \\n 1  MAXTREES     2000  output"
        },
        "output":  {
            "tpl": "//Parameters for the coalescence simulation program : fsimcoal2.exe\\n2 samples to simulate :\\n//Population effective sizes (number of genes)\\n10000\\n8000\\n//Samples sizes and samples age \\n5\\n5\\n//Growth rates: negative growth implies population expansion\\n0\\n0\\n//Number of migration matrices: 0 implies no migration between demes\\n0\\n//Historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index\\n1 historical event\\n2000 1 0 1 7000 0 0\\n//Number of independent loci [chromosome]\\n1 0\\n//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci\\n1\\n//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters\\nFREQ 1 0 2.5e-8 OUTEXP",
            "est": "// Search ranges and rules file\\n// ****************************\\n\\n[PARAMETERS]\\n//#isInt? #name   #dist.#min  #max \\n1  NPOP1       logunif  10   1e7   output\\n1  NPOP2       logunif  10   1e7   output\\n1  NANC        logunif  10   1e7   output \\n1  TDIV        unif     10   1e5   output \\n[RULES]\\n[COMPLEX PARAMETERS]\\n0  RESIZE0   = NANC/NPOP1      hide",
            "par": "[PARAMETERS]\\n1  REPLICATES   100   output\\n1  MAXTREES     2000  output"
        }
    },
        {
        "input": {
            "goal": "analyze population dynamics across three populations with divergence",
            "params": {"NPOP1": 3000, "NPOP2": 5000, "NPOP3": 4000, "TDIV1": 800, "TDIV2": 600},
            "task": "Generate .tpl, .est, and .par files",
            "par_template": "[PARAMETERS] \\n1  REPLICATES   100   output \\n 1  MAXTREES     1500  output"
        },
        "output": {
            "tpl": "//Parameters for the coalescence simulation program : fsimcoal2.exe\\n3 samples to simulate :\\n//Population effective sizes (number of genes)\\n3000\\n5000\\n4000\\n//Samples sizes and samples age \\n5\\n5\\n5\\n//Growth rates: negative growth implies population expansion\\n0\\n0\\n0\\n//Number of migration matrices: 0 implies no migration between demes\\n0\\n//Historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index\\n2 historical events\\n800 1 2 1 4500 0 0\\n600 2 3 1 4000 0 0\\n//Number of independent loci [chromosome]\\n1 0\\n//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci\\n1\\n//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters\\nFREQ 1 0 2.5e-8 OUTEXP",
            "est": "// Search ranges and rules file\\n// ****************************\\n\\n[PARAMETERS]\\n1  NPOP1       logunif  10   1e7   output\\n1  NPOP2       logunif  10   1e7   output\\n1  NPOP3       logunif  10   1e7   output\\n1  TDIV1       unif     10   1e5   output\\n1  TDIV2       unif     10   1e5   output\\n[RULES]\\n[COMPLEX PARAMETERS]\\n0  RESIZE0   = (NPOP1 + NPOP2) / 2      hide",
            "par": "[PARAMETERS]\\n1  REPLICATES   100   output\\n1  MAXTREES     1500  output"
        }
    },
    {
        "input": {
            "goal": "explore dynamics across four populations",
            "params": {"NPOP1": 6000, "NPOP2": 5000, "NPOP3": 8000, "NPOP4": 7000, "TDIV": 1200},
            "task": "Generate .tpl, .est, and .par files",
            "par_template": "[PARAMETERS] \\n 1  REPLICATES   200   output \\n 1  MAXTREES     3000  output"
        },
        "output": {
            "tpl": "//Parameters for the coalescence simulation program : fsimcoal2.exe\\n4 samples to simulate :\\n//Population effective sizes (number of genes)\\n6000\\n5000\\n8000\\n7000\\n//Samples sizes and samples age \\n5\\n5\\n5\\n5\\n//Growth rates: negative growth implies population expansion\\n0\\n0\\n0\\n0\\n//Number of migration matrices: 0 implies no migration between demes\\n0\\n//Historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index\\n1 historical event\\n1200 1 4 1 7500 0 0\\n//Number of independent loci [chromosome]\\n1 0\\n//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci\\n1\\n//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters\\nFREQ 1 0 2.5e-8 OUTEXP",
            "est": "// Search ranges and rules file\\n// ****************************\\n\\n[PARAMETERS]\\n1  NPOP1       logunif  10   1e7   output\\n1  NPOP2       logunif  10   1e7   output\\n1  NPOP3       logunif  10   1e7   output\\n1  NPOP4       logunif  10   1e7   output\\n1  TDIV        unif     10   1e5   output\\n[RULES]\\n[COMPLEX PARAMETERS]\\n0  RESIZE0   = (NPOP1 + NPOP2 + NPOP3) / 3      hide",
            "par": "[PARAMETERS]\\n1  REPLICATES   200   output\\n1  MAXTREES     3000  output"
        }
    },
    
]




# 权限控制，登录页面
# RAG 加了直接扔进embbeding
# {
#         "input": {
#             "goal": "simulate a large-scale demographic model across 26 populations",
#             "params": {"NPOPi": i * 1000 + 500 for i in range(1, 27)},
#             "task": "Generate .tpl, .est, and .par files",
#             "par_template": """
#             [PARAMETERS]
#             1  REPLICATES   500   output
#             1  MAXTREES     10000  output
#             """
#         },
#         "output": {
#             "tpl": """
#             //Parameters for the coalescence simulation program : fsimcoal2.exe
#             26 samples to simulate :
#             //Population effective sizes (number of genes)
#             """ + "\n".join([str(i * 1000 + 500) for i in range(1, 27)]) + """
#             //Samples sizes and samples age 
#             """ + "\n".join(["5"] * 26) + """
#             //Growth rates: negative growth implies population expansion
#             """ + "\n".join(["0"] * 26) + """
#             //Number of migration matrices: 0 implies no migration between demes
#             0
#             //Historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
#             1 historical event
#             3000 1 26 1 8000 0 0
#             //Number of independent loci [chromosome]
#             1 0
#             //Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
#             1
#             //per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
#             FREQ 1 0 2.5e-8 OUTEXP
#             """,
#             "est": """
#             // Search ranges and rules file
#             // ****************************

#             [PARAMETERS]
#             """ + "\n".join([f"1  NPOP{i}       logunif  10   1e7   output" for i in range(1, 27)]) + """
#             1  TDIV        unif     10   1e5   output

#             [RULES]

#             [COMPLEX PARAMETERS]
#             0  RESIZE0   = SUM(NPOPi)/26      hide
#             """,
#             "par": """
#             [PARAMETERS]
#             1  REPLICATES   500   output
#             1  MAXTREES     10000  output
#             """
#         }
#     }