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

   原版VIC模型采用的基于温度阈值的雨雪划分方法会严重低估降雪量。大量研究表明，综合考虑相对湿度、海拔和气温等参数可以更准确地预测降水相态。为此，我们对VIC模型进行了改进，新增了基于Ding et al. (2014)的参数化方案。该方案通过整合每个网格单元的相对湿度、海拔和湿球温度等地理和气象信息，能够更精确地进行雨雪划分。<br>
   修改后的模型同时保留了原版的温度阈值方案和新的Ding et al. (2014)参数化方案，用户可在运行模型前通过全局配置文件选择所需的雨雪划分方案。<br>
   
   （注：本次修改还调整了模型处理气象强迫数据的计算逻辑。原版VIC在输入降雨和降雪数据但缺少降水数据时，会将降雨和降雪相加作为降水总量，然后调用calc_rainonly函数重新计算降雨和降雪量。改进后的模型则支持直接基于输入的降雨和降雪数据进行水文计算：模型会首先检查输入的气候强迫数据是否包含降雨和降雪数据，若包含则直接使用这些数据进行后续计算；若仅有降水数据，则调用calc_rainonly函数计算降雨和降雪量。需要注意的是，由于Ding et al.的参数化方案是基于中国区域的日尺度数据开发的，在处理非中国区域或者亚日尺度数据时可能存在一定不确定性。因此，我们建议用户直接输入降雨和降雪数据进行计算。*CMFD-2.0*数据集包含VIC-5模型运行所需的7个基本变量，同时提供了降雨和降雪数据，是推荐的输入数据源。）
   
   下面的字段需要在全局文件中指定：TEMP_TH_TYPE
   
   1. 原版基于温度阈值的VIC需要MAX_SNOW_TEMP和MIN_RAIN_TEMP两个参数（默认分别为0.5和-0.5），这两个参数目前已经修改为在土壤文件中（具体修改见read_soilparam.c）。<br>
   ```
   TEMP_TH_TYPE VIC_CLASSIC
   ```
   2. 基于Ding et al. (2014)的参数化方案（推荐，默认设置）<br> 
   ```
   TEMP_TH_TYPE VIC_DING
   ```
3. 水文响应单元<br>

   参考了VIC-GL模型中的相关概念，但对energy_bal_struct，snow_data_struct等数组的结构进行了修改。
   
## References
   [1]	Ding B, Yang K, Qin J, et al. The dependence of precipitation types on surface elevation and meteorological conditions and its parameterization[J]. Journal of hydrology, 2014, 513: 154-163.



   


