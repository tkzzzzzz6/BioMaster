from agents.Biomaster import Biomaster
from langchain_core.messages import HumanMessage
import yaml
import json
import os
import sys

# 从YAML配置文件加载配置
def load_config(config_path='config.yaml'):
    """加载YAML配置文件"""
    with open(config_path, 'r', encoding='utf-8') as file:
        return yaml.safe_load(file)

# 主函数
if __name__ == "__main__":
    # 解析命令行参数，获取配置文件路径
    config_path = sys.argv[1] if len(sys.argv) > 1 else 'config.yaml'
    
    # 如果配置文件不存在，报错并退出
    if not os.path.exists(config_path):
        print(f"错误: 配置文件 '{config_path}' 不存在")
        sys.exit(1)
    
    # 加载配置
    config = load_config(config_path)
    
    # 从配置中获取API设置
    api_key = config['api']['main']['key']
    base_url = config['api']['main']['base_url']
    embedding_api_key = config['api']['embedding']['key']
    embedding_base_url = config['api']['embedding']['base_url']
    
    # 从配置中获取模型设置
    Model = config['models']['main']
    tool_model = config['models']['tool']
    embedding_model = config['models']['embedding']
    
    # 从配置中获取Biomaster设置
    excutor = config['biomaster']['executor']
    ids = config['biomaster']['id']
    generate_plan = config['biomaster'].get('generate_plan', True)  # 默认为True
    
    # 从配置中获取数据和目标
    datalist = config['data']['files']
    goal = config['data']['goal']

    # 初始化Biomaster实例
    manager = Biomaster(
        api_key, 
        base_url,
        excutor=excutor,
        id=ids,
        Model=Model,
        embedding_model=embedding_model,
        tool_model=tool_model,
        embedding_base_url=embedding_base_url,
        embedding_api_key=embedding_api_key
    )

    # 根据配置决定是否执行PLAN
    if generate_plan:
        manager.execute_PLAN(goal, datalist)
        print("**********************************************************")

    # 执行TASK
    PLAN_results_dict = manager.execute_TASK(datalist)
    print(PLAN_results_dict)
