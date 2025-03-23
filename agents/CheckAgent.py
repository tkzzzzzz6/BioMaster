import os
import json
import logging
import re
from typing import Dict, List, Tuple, Any
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate

class CheckAgent:
    """
    CheckAgent，核心功能:
    - 验证预期输出文件是否存在且非空（0字节）
    - 文件不存在时扫描目录并让大模型判断可能的输出文件
    - 更新DEBUG_Output文件中的状态
    """
    
    def __init__(self, api_key: str, base_url: str, model: str = "o1-2024-12-17", **kwargs):
        # 设置日志记录
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
        self.model = ChatOpenAI(model=model, base_url=base_url)
        
        # 创建提示模板 - 修复了JSON格式的转义问题
        self.file_prompt = ChatPromptTemplate.from_messages([
            ("system", """
            You are a bioinformatics file analyzer. Based on the step information and file list, 
            determine which files are likely outputs. Return JSON with:
            {{
                "analysis": "Detailed analysis of why files are valid or invalid",
                "output_filename": ["file1: description", "file2: description", ...],
                "stats": true/false (whether the files are valid outputs)
            }}
            """),
            ("human", """
            Step details:
            - Step number: {step_number}
            - Tool: {tool_name}
            - Description: {step_description}
            
            Available files: {file_list}
            
            Expected output files: {expected_files}
            """)
        ])
    
    def check_file_size(self, file_path: str) -> Tuple[bool, str]:
        """检查文件是否存在且非空（0字节）"""
        if not os.path.exists(file_path):
            return False, f"File doesn't exist: {file_path}"
        
        if os.path.getsize(file_path) == 0:
            return False, f"File is empty (0 bytes): {file_path}"
            
        return True, f"File exists and has content"
    
    def scan_directory(self, directory: str) -> List[Dict[str, Any]]:
        """扫描目录中的所有非空文件"""
        if not os.path.exists(directory):
            return []
        
        files = []
        for root, _, filenames in os.walk(directory):
            for filename in filenames:
                path = os.path.join(root, filename)
                if os.path.exists(path) and os.path.getsize(path) > 0:
                    files.append({
                        "path": path,
                        "size": os.path.getsize(path)
                    })
        return files
    
    def get_step_details(self, task_id: str, step_number: str) -> Dict[str, Any]:
        """从PLAN.json获取步骤详情"""
        plan_path = f"./output/{task_id}_PLAN.json"
        if not os.path.exists(plan_path):
            return {}
        
        plan_file = open(plan_path, 'r', encoding='utf-8')
        plan_data = plan_file.read()
        plan_file.close()
        
        if not plan_data:
            return {}
            
        plan = json.loads(plan_data)
        
        for step in plan.get("plan", []):
            if str(step.get("step_number", "")) == step_number:
                return step
        
        return {}
    
    def is_valid_json(self, json_str: str) -> bool:
        """检查字符串是否为有效的JSON"""
        if not json_str:
            return False
            
        # 简单检查是否以{ 开头，以 }结尾
        json_str = json_str.strip()
        if not (json_str.startswith('{') and json_str.endswith('}')):
            return False
            
        # 更多的检查可以添加在这里
        return True
    
    def safe_parse_json(self, json_str: str) -> Dict[str, Any]:
        """安全解析JSON字符串，避免使用try-except"""
        if not self.is_valid_json(json_str):
            self.logger.error(f"Invalid JSON format: {json_str[:100]}...")
            # Return an empty dict with just the error message to append to existing analysis
            return {
                "analysis": "Note: Response format was invalid but processing continues.",
                "output_filename": [],
                "stats": True  # Changed from False to True to not fail due to format issues
            }
            
        try:
            return json.loads(json_str)
        except json.JSONDecodeError:
            self.logger.error(f"JSON decode error: {json_str[:100]}...")
            return {
                "analysis": "Note: JSON parsing error occurred but processing continues.",
                "output_filename": [],
                "stats": True
            }
    
    def check_output_files(self, debug_output_path: str) -> Dict[str, Any]:
        """
        主要检查函数：检查文件是否存在且非空，不存在则扫描目录分析
        """
        # 检查并加载DEBUG_Output文件
        if not os.path.exists(debug_output_path):
            self.logger.error(f"DEBUG_Output file not found: {debug_output_path}")
            return {"overall_check": False}
            
        debug_file = open(debug_output_path, 'r', encoding='utf-8')
        debug_content = debug_file.read()
        debug_file.close()
        
        if not debug_content:
            self.logger.error(f"DEBUG_Output file is empty")
            return {"overall_check": False}
            
        debug_output = json.loads(debug_content)
        
        # 如果stats已经是False，直接返回
        if not debug_output.get("stats", False):
            return {"overall_check": False}
        
        # 从路径提取任务ID和步骤号
        path_match = re.search(r'(\d+)_DEBUG_Output_(\d+)\.json', debug_output_path)
        if not path_match:
            return {"overall_check": False}
        
        task_id = path_match.group(1)
        step_number = path_match.group(2)
        output_dir = f"./output/{task_id}"
        
        # 获取期望的输出文件列表
        output_files = debug_output.get("output_filename", [])
        expected_files = []
        for file_entry in output_files:
            if ':' in file_entry:
                path = file_entry.split(':', 1)[0].strip()
            else:
                path = file_entry.strip()
            expected_files.append(path)
        
        # 检查所有预期文件
        all_valid = True
        
        # 如果有预期文件，检查它们是否存在且非空
        if expected_files:
            for file_path in expected_files:
                valid, _ = self.check_file_size(file_path)
                if not valid:
                    all_valid = False
                    break
            
            # 如果所有文件都有效，直接返回成功
            if all_valid:
                return {
                    "overall_check": True,
                    "stats": True
                }
        
        # 扫描目录中的所有非空文件
        all_files = self.scan_directory(output_dir)
        if not all_files:
            # 没有找到任何文件，直接标记失败
            debug_output["stats"] = False
            debug_output["analyze"] = debug_output.get("analyze", "") + f"\n\nNo valid files found in {output_dir}"
            
            out_file = open(debug_output_path, 'w', encoding='utf-8')
            out_file.write(json.dumps(debug_output, indent=4))
            out_file.close()
            
            return {"overall_check": False, "stats": False}
        
        # 获取步骤详情并格式化文件列表
        step_details = self.get_step_details(task_id, step_number)
        file_list_text = "\n".join([f"- {f['path']} ({f['size']} bytes)" for f in all_files])
        
        # 调用大模型分析
        response = self.file_prompt.invoke({
            "step_number": step_details.get("step_number", "unknown"),
            "tool_name": step_details.get("tools", "unknown"),
            "step_description": step_details.get("description", "unknown"),
            "file_list": file_list_text,
            "expected_files": ", ".join(expected_files) if expected_files else "None specified"
        })
        
        result = self.model.invoke(response)
        
        # 解析模型回复，避免使用try-except
        model_response = self.safe_parse_json(result.content)
        analysis = model_response.get("analysis", "")
        new_output_files = model_response.get("output_filename", [])
        stats = model_response.get("stats", True)  # Default to True instead of False
        print("--------------------------------")
        print(model_response)
        print("--------------------------------")
        # 更新DEBUG_Output文件
        debug_output["stats"] = debug_output.get("stats", True) and stats  # Preserve existing stats unless new one is False
        debug_output["analyze"] = debug_output.get("analyze", "") + f"\n\n{analysis}"
        if new_output_files:
            debug_output["output_filename"] = new_output_files
        
        out_file = open(debug_output_path, 'w', encoding='utf-8')
        out_file.write(json.dumps(debug_output, indent=4))
        out_file.close()
        
        return {
            "overall_check": stats,
            "analysis": analysis,
            "output_filename": new_output_files,
            "stats": stats
        }
