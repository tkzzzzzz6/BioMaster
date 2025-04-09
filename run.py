from agents.Biomaster import Biomaster
from langchain_core.messages import HumanMessage
import json
# Example of using the agent
if __name__ == "__main__":

    config = {"configurable": {"thread_id": "abc124"}}
    api_key = 'sk-YCGdylUyNQjXPRB246CbEb78BcA24578B821728f2d8358Ba'
    base_url = 'https://api.bltcy.ai/v1'
    # base_url = 'https://api.siliconflow.cn/v1'
    # you can choose other base url
    manager = Biomaster(api_key, base_url,excutor=True,id='004')
    
    # datalist=["./data/1000GP_pruned.bed: SNP file in bed format",
    #             "./data/1000GP_pruned.bim: snp info associate with the bed format",
    #             "./data/1000GP_pruned.fam: data information, the first col is population, the second is sample ID",]
    # goal='please help me do ROH.'
    # datalist=['./data/SRR620204.fastq.gz: chip-seq data for Ring1B',
    #             './data/SRR620205.fastq.gz: chip-seq data for cbx7',
    #             './data/SRR620206.fastq.gz: chip-seq data for SUZ12',
    #             './data/SRR620208.fastq.gz: chip-seq data for IgGold',
    #             './data/mm39.ncbiRefSeq.gtf: genome annotations mouse',
    #             './data/mm39.fa: mouse genome',]
    # goal='call peaks for protein cbx7 with IgGold as control'
    datalist=['./data/SRR620204.fastq.gz: chip-seq data for Ring1B',
                './data/SRR620205.fastq.gz: chip-seq data for cbx7',
                './data/SRR620206.fastq.gz: chip-seq data for SUZ12',
                './data/SRR620208.fastq.gz: chip-seq data for IgGold',
                './data/mm39.ncbiRefSeq.gtf: genome annotations mouse',
                './data/mm39.fa: mouse genome',]
    goal='call peaks for protein cbx7 with IgGold as control'
    manager.execute_PLAN(goal,datalist)
    print("**********************************************************")

    PLAN_results_dict = manager.execute_TASK(datalist)
    print(PLAN_results_dict)



