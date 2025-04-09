from langchain_openai import ChatOpenAI
from langgraph.checkpoint.sqlite import SqliteSaver
from langchain_core.prompts import ChatPromptTemplate
import bs4
from langchain_community.document_loaders import WebBaseLoader
from langchain.prompts import ChatPromptTemplate, FewShotChatMessagePromptTemplate
from langchain_community.document_loaders.text import TextLoader
from langchain_openai import ChatOpenAI
from langchain_core.output_parsers import StrOutputParser
from langchain import hub
from langchain.agents import AgentExecutor, create_tool_calling_agent
from langchain_community.tools import ShellTool
import os
from openai import OpenAI

def Json_Format_Agent(text, api_key, base_url, model="o3-mini", return_tokens=False):
    # 直接使用 OpenAI 客户端而不是 LangChain
    client = OpenAI(api_key=api_key, base_url=base_url)
    
    # 创建一个用于修正 JSON 的提示模板
    messages = [
        {
            "role": "system",
            "content": """You are an expert in JSON formatting.
            You cannot reply to anything other than json format.
            Do not reply to anything outside of json.
            You don't need to put json in the code block.
            You must strictly observe all the above rules.
            Do not delete any extra Spaces."""
        },
        {
            "role": "user",
            "content": f"Given the following potentially malformed JSON string, provide a correctly formatted JSON string:\n\n{text}\n\nCorrect JSON:"
        }
    ]
    
    # 调用API
    response = client.chat.completions.create(
        model=model,
        messages=messages
    )
    
    formatted_text = response.choices[0].message.content
    
    if return_tokens:
        input_tokens = response.usage.prompt_tokens
        output_tokens = response.usage.completion_tokens
        return formatted_text, input_tokens, output_tokens
    
    return formatted_text