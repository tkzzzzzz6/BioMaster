from agents.Biomaster import Biomaster
from langchain_core.messages import HumanMessage
import json
# Example of using the agent
if __name__ == "__main__":

    config = {"configurable": {"thread_id": "abc124"}}

    # support main model:
    # o3-mini, o1, gpt-4o, o3-mini-2025-01-31, o1-2024-12-17
    # claude-3-7-sonnet-thinking, claude-3-7-sonnet-20250219, claude-3-5-sonnet-20241022
    # DeepSeek-V3, Deepseek-R1
    # Qwen/QWQ-32B 
    # LLAMA3-70B
    # All other models can be tried, but it is suggested that if the main model chooses a better model, 
    # the best model tested at present is o3-mini-2025-01-31

    # support tool model:
    # All LLMS can be tried, and you can choose some small models here.

    # support emmbedding model
    # text-embedding-004,
    # text-embedding-3-large, text-embedding-3-small,
    # text-embedding-ada-002
    # BAAI/bge-m3

    # suggest base url:
    # https://api.bltcy.ai/v1
    # https://gpt-api.hkust-gz.edu.cn/v1
    # https://dashscope.aliyuncs.com/compatible-mode/v1
    # https://api.openai.com/v1
    # https://api.siliconflow.cn/v1
    # https://sg.uiuiapi.com/v1

    api_key = ''
    base_url = 'https://api.siliconflow.cn/v1/'
    embedding_api_key = ''
    embedding_base_url = 'https://api.siliconflow.cn/v1/'
    embedding_model='BAAI/bge-m3'
    Model=tool_model = "deepseek-ai/DeepSeek-V3"
    tool_model = "deepseek-ai/DeepSeek-V3"
    # you can choose other base url
    manager = Biomaster(api_key, base_url,excutor=True,id='001',Model=Model,embedding_model=embedding_model,tool_model=tool_model,embedding_base_url=embedding_base_url,embedding_api_key=embedding_api_key)

    datalist=[ './data/rnaseq_1.fastq.gz: RNA-Seq read 1 data (left read)',
            './data/rnaseq_2.fastq.gz: RNA-Seq read 2 data (right read)',
            './data/minigenome.fa: small genome sequence consisting of ~750 genes.',]
    goal='please do WGS/WES data analysis Somatic SNV+indel calling.'
    manager.execute_PLAN(goal,datalist)
    print("**********************************************************")

    PLAN_results_dict = manager.execute_TASK(datalist)
    print(PLAN_results_dict)
