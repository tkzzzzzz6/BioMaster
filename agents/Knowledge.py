import os
import json
import tempfile
import pytesseract
import pdfplumber
import requests
import io
import re
from PIL import Image
from pdf2image import convert_from_path
from bs4 import BeautifulSoup
from urllib.parse import urlparse, urljoin

from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

# 引入你自定义的 JSON 修复函数
from .ToolAgent import Json_Format_Agent

def extract_text_from_pdf(pdf_path: str) -> str:
    """
    先尝试用 pdfplumber 提取文本；如果没有文本（扫描版），
    则使用 pdf2image + pytesseract 来进行 OCR。
    """
    extracted_text_pages = []
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            page_text = page.extract_text() or ""
            extracted_text_pages.append(page_text)
    all_text = "\n".join(extracted_text_pages).strip()

    if all_text:
        return all_text

    # 如果 pdfplumber 没取到文本，就 OCR
    ocr_text_pages = []
    with tempfile.TemporaryDirectory() as tempdir:
        pages = convert_from_path(pdf_path, dpi=300, output_folder=tempdir)
        for page_image in pages:
            text = pytesseract.image_to_string(page_image)
            ocr_text_pages.append(text)
    return "\n".join(ocr_text_pages)

def extract_text_from_url(url: str, process_images: bool = True, api_key: str = None) -> str:
    """
    从网页URL提取文本内容，可选择性处理图片
    
    参数:
        url: 网页地址
        process_images: 是否处理网页中的图片
        api_key: 用于视觉模型的API密钥
    """
    try:
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()  # 确保请求成功
        
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # 移除script和style元素
        for script in soup(["script", "style"]):
            script.decompose()
            
        # 获取文本
        text = soup.get_text(separator='\n')
        
        # 清理文本
        lines = (line.strip() for line in text.splitlines())
        chunks = (phrase.strip() for line in lines for phrase in line.split("  "))
        text = '\n'.join(chunk for chunk in chunks if chunk)
        
        # 如果启用了图片处理且提供了API密钥
        if process_images and api_key:
            image_texts = extract_text_from_images(soup, url, api_key)
            if image_texts:
                text += "\n\n--- 图片内容 ---\n" + "\n\n".join(image_texts)
        
        return text
    except Exception as e:
        return f"Error extracting text from URL: {str(e)}"

def extract_text_from_images(soup, base_url, api_key):
    """提取网页中的图片并进行OCR或视觉分析"""
    image_texts = []
    
    # 查找所有图片标签
    for img in soup.find_all('img'):
        try:
            # 获取图片URL
            img_url = img.get('src')
            if not img_url:
                continue
                
            # 将相对URL转为绝对URL
            img_url = urljoin(base_url, img_url)
            
            # 判断是否是真实图片 (排除图标等)
            if img.get('width') and int(img.get('width')) < 100:
                continue
                
            # 下载图片
            img_response = requests.get(img_url, stream=True)
            if img_response.status_code != 200:
                continue
                
            # 打开图片
            img_content = Image.open(io.BytesIO(img_response.content))
            
            # 进行OCR处理
            text = pytesseract.image_to_string(img_content)
            
            # 只保留有实质内容的OCR结果
            if len(text.strip()) > 10:  # 假设少于10个字符的是噪音
                image_texts.append(text)
                
            # 可选: 也可以替换/增加为视觉模型处理
            # vision_text = analyze_image_with_vision_model(img_content, api_key)
            # if vision_text:
            #     image_texts.append(vision_text)
            
        except Exception as e:
            continue  # 忽略处理单张图片的错误，继续处理其他图片
            
    return image_texts

def analyze_image_with_vision_model(image, api_key):
    """
    使用视觉语言模型分析图片内容
    这是一个示例实现，需要根据实际使用的API进行调整
    """
    try:
        # 这里是集成视觉语言模型API的示例代码
        # 例如可以使用OpenAI的GPT-4 Vision或其他视觉语言模型API
        
        # 示例: 如果使用OpenAI
        # 首先需要转换图像为base64格式
        import base64
        from io import BytesIO
        
        buffered = BytesIO()
        image.save(buffered, format="JPEG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        # 然后调用API
        # response = openai.chat.completions.create(
        #     model="gpt-4-vision-preview",
        #     messages=[
        #         {
        #             "role": "user",
        #             "content": [
        #                 {"type": "text", "text": "What's in this image? Extract any text visible and describe key elements."},
        #                 {"type": "image_url", "image_url": {"url": f"data:image/jpeg;base64,{img_str}"}}
        #             ]
        #         }
        #     ],
        #     max_tokens=300
        # )
        # return response.choices[0].message.content
        
        # 由于这部分需要实际的API集成，返回占位符
        return "视觉语言模型功能待集成"
        
    except Exception as e:
        return f"Error analyzing image: {str(e)}"

def summarize_content_as_workflow(content_path: str, api_key: str, base_url: str, process_images: bool = True) -> str:
    """
    1. 从PDF或网页中提取文本
    2. 使用GPT对文本进行总结，输出指定JSON结构
    3. 调用Json_Format_Agent修复JSON确保合法
    
    参数:
        content_path: PDF文件路径或网页URL
        api_key: OpenAI API密钥
        base_url: OpenAI API基础URL
        process_images: 处理网页图片（仅对URL有效）
    """
    # 判断输入是PDF文件还是URL
    is_url = bool(urlparse(content_path).scheme)
    
    if is_url:
        all_text = extract_text_from_url(content_path, process_images, api_key)
    else:
        all_text = extract_text_from_pdf(content_path)

    template_for_summary = """
    You are an expert bioinformatics workflow specialist. Your task is to read the provided scientific paper and extract the methodological workflows and techniques used in the research.
    
    Please analyze the paper and create a detailed step-by-step workflow with the following requirements:
    
    1. Number each step sequentially (Step 1, Step 2, Step 3, etc.)
    2. Break down the workflow into fine-grained steps, with each step using preferably only ONE tool
    3. For each step, clearly specify:
       - Description: Brief explanation of what this step accomplishes
       - Tool Used: The specific software/tool used in this step
       - Input Required: What data or files are needed as input
       - Expected Output: What results or files are produced
    4. Be as comprehensive as possible, capturing all methodological details
    
    Format your response as a valid JSON according to the following structure:
    
    {{
        "content": "This is a workflow for [Brief description of the workflow]\n\nStep 1: [Step Name]\nDescription: [Brief description of the step]\nInput Required: [Input data/files needed]\nExpected Output: [Output data/files produced]\nTool Used: [Specific tool name].\n\nStep 2: [Step Name]\nDescription: [Brief description of the step]\nInput Required: [Input data/files needed]\nExpected Output: [Output data/files produced]\nTool Used: [Specific tool name].\n\n[Continue with remaining steps...]",
        "metadata": {{
            "source": "workflow",
            "page": 1
        }}
    }}
    
    Paper text to analyze:
    {text}
    """

    os.environ['OPENAI_API_KEY'] = api_key

    summary_prompt = ChatPromptTemplate.from_template(template_for_summary)
    llm = ChatOpenAI(model="o3-mini-2025-01-31", base_url=base_url)
    parser = StrOutputParser()
    chain = summary_prompt | llm | parser

    # 调用模型得到初步 JSON
    raw_json_output = chain.invoke({"text": all_text})
    # 修复 JSON
    corrected_json = Json_Format_Agent(raw_json_output, api_key=api_key, base_url=base_url)
    return corrected_json

# 保留原函数名作为别名，保持向后兼容
summarize_pdf_as_workflow = summarize_content_as_workflow

def append_summary_to_plan_knowledge(corrected_json: str, plan_knowledge_path: str = "./doc/Plan_Knowledge.json"):
    # 如果文件不存在，先创建一个空的 JSON 数组文件
    if not os.path.exists(plan_knowledge_path):
        with open(plan_knowledge_path, 'w', encoding='utf-8') as f:
            f.write("[]")

    # 读取原有 JSON
    with open(plan_knowledge_path, 'r', encoding='utf-8') as f:
        try:
            existing_data = json.load(f)
        except json.JSONDecodeError:
            existing_data = []

    if not isinstance(existing_data, list):
        existing_data = [existing_data]

    # 将字符串形式的 corrected_json 转成 Python 对象
    try:
        new_summary_obj = json.loads(corrected_json)
    except json.JSONDecodeError:
        print("传入的 corrected_json 不是合法 JSON 字符串，请检查。")
        return

    # 追加到列表
    existing_data.append(new_summary_obj)

    # 写回文件
    with open(plan_knowledge_path, 'w', encoding='utf-8') as f:
        json.dump(existing_data, f, ensure_ascii=False, indent=4)

def generate_workflow_markdown(content_path: str, api_key: str, base_url: str, output_path: str = "plan_knowledge.md", process_images: bool = True) -> str:
    """
    从文档内容生成Markdown格式的工作流描述，并保存到文件
    
    参数:
        content_path: PDF文件路径或网页URL
        api_key: OpenAI API密钥
        base_url: OpenAI API基础URL
        output_path: 输出的Markdown文件路径
        process_images: 处理网页图片（仅对URL有效）
    
    返回:
        生成的Markdown内容字符串
    """
    # 判断输入是PDF文件还是URL
    is_url = bool(urlparse(content_path).scheme)
    
    if is_url:
        all_text = extract_text_from_url(content_path, process_images, api_key)
    else:
        all_text = extract_text_from_pdf(content_path)

    template_for_summary = """
    You are an expert bioinformatics workflow specialist. Your task is to read the provided scientific paper and extract the methodological workflows and techniques used in the research.
    
    Please analyze the paper and create a detailed step-by-step workflow with the following requirements:
    
    1. Number each step sequentially (Step 1, Step 2, Step 3, etc.)
    2. Break down the workflow into fine-grained steps, with each step using preferably only ONE tool
    3. For each step, clearly specify:
       - Description: Brief explanation of what this step accomplishes
       - Tool Used: The specific software/tool used in this step
       - Input Required: What data or files are needed as input
       - Expected Output: What results or files are produced
    4. Be as comprehensive as possible, capturing all methodological details
    
    Format your response in Markdown format starting with:
    
    # Workflow for [Brief description of the workflow]
    
    ## Step 1: [Step Name]
    **Description:** [Brief description of the step]
    **Input Required:** [Input data/files needed]
    **Expected Output:** [Output data/files produced]
    **Tool Used:** [Specific tool name]
    
    ## Step 2: [Step Name]
    **Description:** [Brief description of the step]
    **Input Required:** [Input data/files needed]
    **Expected Output:** [Output data/files produced]
    **Tool Used:** [Specific tool name]
    
    [Continue with remaining steps...]
    
    Paper text to analyze:
    {text}
    """

    os.environ['OPENAI_API_KEY'] = api_key

    summary_prompt = ChatPromptTemplate.from_template(template_for_summary)
    llm = ChatOpenAI(model="o3-mini-2025-01-31", base_url=base_url)
    parser = StrOutputParser()
    chain = summary_prompt | llm | parser

    # 调用模型生成Markdown
    markdown_output = chain.invoke({"text": all_text})
    
    # 保存Markdown到文件
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(markdown_output)
        
    print(f"已生成工作流Markdown文件: {output_path}")
    return markdown_output

def convert_markdown_to_json(markdown_path: str, json_path: str = "./doc/Plan_Knowledge.json"):
    """
    将Markdown格式的工作流转换为JSON格式并保存
    
    参数:
        markdown_path: Markdown文件路径
        json_path: 输出的JSON文件路径
    """
    try:
        # 读取Markdown文件
        with open(markdown_path, 'r', encoding='utf-8') as f:
            markdown_content = f.read()
        
        # 提取工作流标题
        title_match = re.search(r'# Workflow for (.*?)(?:\n|$)', markdown_content)
        workflow_title = title_match.group(1) if title_match else "Unnamed Workflow"
        
        # 构建文本内容，转换Markdown为纯文本格式
        text_content = f"This is a workflow for {workflow_title}\n\n"
        
        # 提取每个步骤
        step_pattern = r'## Step (\d+): (.*?)(?:\n|$)(.*?)(?=## Step \d+:|$)'
        steps = re.findall(step_pattern, markdown_content, re.DOTALL)
        
        for step_num, step_name, step_content in steps:
            # 提取步骤详情
            description = re.search(r'\*\*Description:\*\* (.*?)(?:\n|$)', step_content)
            input_req = re.search(r'\*\*Input Required:\*\* (.*?)(?:\n|$)', step_content)
            output = re.search(r'\*\*Expected Output:\*\* (.*?)(?:\n|$)', step_content)
            tool = re.search(r'\*\*Tool Used:\*\* (.*?)(?:\n|$)', step_content)
            
            # 构建步骤文本
            step_text = f"Step {step_num}: {step_name}\n"
            if description: step_text += f"Description: {description.group(1)}\n"
            if input_req: step_text += f"Input Required: {input_req.group(1)}\n"
            if output: step_text += f"Expected Output: {output.group(1)}\n"
            if tool: step_text += f"Tool Used: {tool.group(1)}\n"
            
            text_content += step_text + "\n"
        
        # 构建JSON对象
        json_obj = {
            "content": text_content,
            "metadata": {
                "source": "workflow",
                "page": 1
            }
        }
        
        # 准备保存到Plan_Knowledge.json
        if not os.path.exists(os.path.dirname(json_path)) and os.path.dirname(json_path):
            os.makedirs(os.path.dirname(json_path))
            
        if os.path.exists(json_path):
            # 读取现有JSON
            with open(json_path, 'r', encoding='utf-8') as f:
                try:
                    existing_data = json.load(f)
                except json.JSONDecodeError:
                    existing_data = []
        else:
            existing_data = []
            
        if not isinstance(existing_data, list):
            existing_data = [existing_data]
            
        # 添加新的JSON对象
        existing_data.append(json_obj)
        
        # 保存到文件
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(existing_data, f, ensure_ascii=False, indent=4)
            
        print(f"已将Markdown转换为JSON并保存至: {json_path}")
        return json.dumps(json_obj, ensure_ascii=False)
    
    except Exception as e:
        print(f"转换过程中出错: {str(e)}")
        return None

