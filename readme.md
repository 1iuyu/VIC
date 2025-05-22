# Variable Infiltration Capacity (VIC) Model
## 这是一个基于VIC-5.1.0的修改版模型
主要的修改内容：
- 增加了对冰川质量和能量平衡计算的支持
- 增加了新的雨雪划分参数化方案
- 增加了水文响应单元提高对模型空间异质性的描述

## 修改说明：
1. 冰川质量和能量平衡计算<br>

   冰川过程的计算参考了PCIC的VIC-GL模型（VIC Glacier），而后者基于的VIC版本太旧，新版VIC在蒸散发等过程都有相关改进，有必要在新版VIC的基础上增加对于冰川过程的计算。

2. 新的雨雪划分参数化方案<br>

基于Ding et al. (2014)的参数化方案，根据每个单元的相对湿度、海拔和湿球温度等地理和气象信息在每个时间步长中进行雨雪划分![image](https://github.com/user-attachments/assets/fc4e91d2-bedf-49e8-8ed0-47eb4bccb52f)


