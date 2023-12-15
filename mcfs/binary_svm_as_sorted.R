setwd("D:/bilibiliR/29、二分类结局变量svm变量筛选")

# 下载需要的R包 -----------------------------------------------------------------
# 检查 e1071 包是否已经安装
if (!require("e1071", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("e1071")
}

# 检查 caret 包是否已经安装
if (!require("caret", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("caret")
}

# 检查 pROC 包是否已经安装
if (!require("pROC", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("pROC")
}

# 检查 ggplot2 包是否已经安装
if (!require("ggplot2", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggplot2")
}

# 加载所需要的R包 ----------------------------------------------------------------
library(e1071)#支持向量机
library(caret)#集成了上百种分类和回归算法
library(pROC)#ROC曲线可视化
library(ggplot2)#绘图所需


# 载入示例数据 ------------------------------------------------------------------

setwd("C:/Users/a/Desktop/表达矩阵")

library(dplyr)
cutoffhold <- 1.5
threshold <- 0.01
heatCD4data<-read.csv("all_diffCD4.csv", header = TRUE)
heatCD4data <- arrange(heatCD4data,by = padj)
#heatCD4 <- heatCD4data$ext_gene[1:150]
# heatCD4data <- heatCD4data[heatCD4data$log2FoldChange > 2,]
heatCD4 <- heatCD4data$ext_gene[heatCD4data$padj < threshold]

heatCD8data<-read.csv("all_diffCD8.csv", header = TRUE)
heatCD8data <- arrange(heatCD8data,by = padj)
# heatCD8data <- heatCD8data[heatCD8data$log2FoldChange > 2,]
heatCD8 <- heatCD8data$ext_gene[heatCD8data$padj < threshold]

heatmodata<-read.csv("all_diffMo.csv", header = TRUE)
heatmodata <- arrange(heatmodata,by = padj)
# heatmodata <- heatCD4data[heatmodata$log2FoldChange > 2,]
heatmo <- heatmodata$ext_gene[heatmodata$padj < threshold]

heatnkdata<-read.csv("all_diffNK.csv", header = TRUE)
heatnkdata <- arrange(heatnkdata,by = padj)
# heatnkdata <- heatnkdata[heatnkdata$log2FoldChange > 2,]
heatnk <- heatnkdata$ext_gene[heatnkdata$padj < threshold]

cd4_cd8 <- intersect(heatCD4,heatCD8)
mo_nk <- intersect(heatmo,heatnk)
selected_gene <- intersect(cd4_cd8,mo_nk)

cutoff_CD4 <- heatCD4data$ext_gene[abs(heatCD4data$log2FoldChange) > cutoffhold]
cutoff_CD8 <- heatCD8data$ext_gene[abs(heatCD8data$log2FoldChange) > cutoffhold]
cutoff_mo <- heatmodata$ext_gene[abs(heatmodata$log2FoldChange) > cutoffhold]
cutoff_nk <- heatCD4data$ext_gene[abs(heatCD4data$log2FoldChange) > cutoffhold]

cutoff_CD4_CD8 <- intersect(cutoff_CD4,cutoff_CD8)
cutoff_mo_nk <- intersect(cutoff_mo,cutoff_nk)
cutoff_all <- intersect(cutoff_CD4_CD8,cutoff_mo_nk)

# heatCD4data <- heatCD4data[heatCD4data$ext_gene %in% selected_gene, ]
# heatCD8data <- heatCD8data[heatCD8data$ext_gene %in% selected_gene, ]
# heatmodata <- heatmodata[heatmodata$ext_gene %in% selected_gene, ]
# heatnkdata <- heatnkdata[heatnkdata$ext_gene %in% selected_gene, ]

heatCD4data <- heatCD4data[, !(names(heatCD4data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatCD8data <- heatCD8data[, !(names(heatCD8data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatmodata <- heatmodata[, !(names(heatmodata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatnkdata <- heatnkdata[, !(names(heatnkdata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]

heatCD4_CD8data <- merge(heatCD4data,heatCD8data, by = "ext_gene")
heatmo_nkdata <- merge(heatmodata,heatnkdata, by = "ext_gene")
heatalldata <- merge(heatCD4_CD8data,heatmo_nkdata, by = "ext_gene")

heatalldata <- heatalldata[heatalldata$ext_gene %in% selected_gene, ]
########为了使图像对比度更高，去掉log2foldchange小于2的值
heatalldata <- heatalldata[heatalldata$ext_gene %in% cutoff_all, ]



write.csv(heatCD4data,file="D:/bulk_control_AS/svm/cd4_count.csv")
write.csv(heatCD8data,file="D:/bulk_control_AS/svm/cd8_count.csv")
write.csv(heatmodata,file="D:/bulk_control_AS/svm/mo_count.csv")
write.csv(heatnkdata,file="D:/bulk_control_AS/svm/nk_count.csv")







# 定义一个函数来根据列名的模式添加新的行
add_control_or_AS <- function(data) {
  # 获取数据框的列名
  col_names <- colnames(data)
  
  # 创建新的行
  new_row <- numeric(length(col_names))
  
  # 遍历列名
  for (i in 1:length(col_names)) {
    if (grepl("^C\\d+", col_names[i])) {  # 列名匹配C紧接着任意数字
      new_row[i] <- "Control"
    } else {
      new_row[i] <- "AS"
    }
  }
  
  # 将新行添加到数据框末尾
  data[nrow(data) + 1, ] <- new_row
  
  # 更改新行的行名为"label"
  rownames(data)[nrow(data)] <- "label" 
  
  
  return(data)
}

# 调用函数来处理数据框
#new_df <- add_control_or_AS(df)



CD4data <- read.csv("D:/bulk_control_AS/svm/cd4_count.csv",row.names = 2)
CD4data$X <- NULL
new_CD4data <- add_control_or_AS(CD4data)
new_CD4data <- t(new_CD4data)


#简单的数据清洗
mydata <- subset(mydata, select = -time)#删除time这一列
#因变量
mydata$status <- factor(ifelse(mydata$status == "Dead",1,0)) #死亡为1，存活为0
#将二分类变量中的字符转变成数字
mydata$gender <- ifelse(mydata$gender == "female", 0, 1)
mydata$stage <- ifelse(mydata$stage == "early", 0, 1)
mydata$type <- ifelse(mydata$type == "squamous carcinom", 0, 1)
mydata$smoker <- ifelse(mydata$smoker == "nonsmoker", 0, 1)
mydata$alcohol_history <- ifelse(mydata$alcohol_history == "No", 0, 1)
#SVM对数值敏感，二分类变量输入会让结果出现很大波动


mydata <- as.data.frame(new_CD4data) 
# SVM算法变量筛选 ---------------------------------------------------------------
table(mydata$label)
X <- mydata[,1:ncol(mydata)-1] # 定义要选择的变量范围
Y <- as.numeric(as.factor(mydata$label)) # 定义预测变量
# ①设置重复交叉验证
control <- trainControl(method = "repeatedcv",   
                        # 交叉验证方法,这里选择重复的k折交叉验证
                        number = 5,     
                        # k值,将数据分成5份,可以根据样本量选择5-10
                        repeats = 3,     
                        # 重复次数,可选2-5次,这里选择1次
                        search = "random"   
                        # 选择重复交叉验证时,每轮选择 folds 的方式,可选"systematic"连续选择与"random"随机选择,这里选择随机     
                        )
  

# ②使用 SVM-RFE 特征选择方法
set.seed(1)   # 为了保持结果的一致性
svm_rfe <- rfe(X,#输入待筛选变量的表达矩阵
               Y,#输入因变量
               #,sizes参数用于指定我们希望选择的特征数量。
               #它决定了RFE算法会考虑选择哪些特征数量的方案。
               #比如sizes = c(1:5, 10)，则RFE算法会考虑选择1个特征,2个特征,3个特征,4个特征,5个特征和10个特征的特征选择方案。
               sizes = 1:20,
               #选择1-10个特征
               rfeControl = rfeControl(functions = caretFuncs,
                                       #使用caret包中的评价指标
                                       method = "repeatedcv",
                                       #使用重复的留一交叉验证
                                       number = 10,
                                       #推荐为10，10折交叉验证
                                       repeats = 5,
                                       #推荐为5，重复5次
                                       verbose = FALSE),
                                      #不打印详细输出
               method = "svmRadial",
               #使用SVM分类器
               trControl = control,
               #trainControl设置
               #preProc = c("center", "scale")
               #特征预处理,标准化
               )

# ③提取前10个最相关的变量
svm_rfe_ranking <- svm_rfe$variables
head(svm_rfe_ranking)

# Overall:特征的总体重要性评分,越高表示该特征对模型性能贡献越大,重要性越高。
# var:特征的名称。
# Variables:特征数量,显示在选取该数量特征时,各特征的重要性评分。这里显示在选择21个特征的情况下,各特征的重要性。
# Resample:如果进行了交叉验证,显示交叉验证的fold信息。这里显示为Fold1.Rep1,表示第一折交叉验证的结果。


varImp(svm_rfe)#查看变量重要性评分

varImp_dataframe <- data.frame(Gene = row.names(varImp(svm_rfe))[1:20],
                               importance = varImp(svm_rfe)[1:20, 1])
                          

# 删除NA行
varImp_dataframe <- na.omit(varImp_dataframe)

# 绘制柱状图
mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')

ggplot(varImp_dataframe, aes(x = reorder(Gene, -importance), y = importance , fill = Gene)) + 
  # 将Gene映射到x轴,按importance降序排序  
  geom_col() +
  # 指定geom_col()为柱状图
  ggtitle("Hub Genes") +
  # 设置标题  
  theme(panel.border = element_blank(),# 删除边框
        axis.text.x = element_text(size = 12, color = "black"),# 修改x轴和y轴文本大小和颜色
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),# 移动x轴标题,增加距离
        plot.title = element_text(margin = margin(b = 20)),# 移动图标题,增加距离  
        panel.grid.major = element_line(color = "grey", size = 0.2)) +# 增加灰色网格线
  # 修改x轴名称为"Gene" ，修改y轴名称为"Importance"
  xlab("Gene") + ylab("Importance") +
  scale_fill_manual(values = mycolors)# 使用自定义的mycolors向量填充颜色  


# 查看重要性前10的基因
top_10_vars <- svm_rfe_ranking$var[1:10]
top_10_vars

# 提取前10变量的表达矩阵
top10_SVM_data <- mydata[,top_10_vars]

# ④提取最优子集中的变量
X_plot = svm_rfe$results$Variables
Y_plot = svm_rfe$results$RMSE
plot(X_plot, Y_plot,  
     xlab="Variable Number",  
     ylab="RMSE (Cross-Validation)",  
     col="#7DEE44",    # 点的颜色
     pch=16,           # 点的形状, ici选择实心圆点       
     cex=1.5,          # 点的大小
     lwd=2,            # 线的宽度
     type="b",         # 同时绘制点与线
     ylim=c(0.54, 0.62)) # y轴范围
#lines(X_plot, Y_plot, col="#DF294C", lwd=2)   # 额外绘制一条红色粗线

abline(h=min(Y_plot), col="skyblue")   # 添加一条水平参考线
grid(col="grey",lwd=1,lty=3)   # 添加网格线 
 
legend("topright",c("Training RMSE","Cross-Validation RMSE"),
       col=c("#7DEE44","#DF294C"),pch=c(16,NA),lwd=2,bg="white")   # 添加图例

# 找到RMSE最小的点
# RMSE是Root Mean Square Error(均方根误差)的缩写,是评价回归模型效果的一个很重要指标。
# RMSE通过测量预测值与实际值的离差来评价模型的准确性,值越小表示模型越准确。
wmin <- which.min(svm_rfe$results$RMSE)
wmin

# 在图上标记RMSE最小的点
points(wmin, svm_rfe$results$RMSE[wmin], col = "orange", pch = 16, cex=2)  
text(wmin, svm_rfe$results$RMSE[wmin],  
     paste0("N=", wmin), pos = 2, col = "orange", cex=2)  

# 提取最优子集中变量的名称
Target_Genes <- svm_rfe$optVariables
Target_Genes

# 提取最优子集的表达矩阵
Best_SVM_data <- mydata[,Target_Genes]



# 支持向量机拟合模型 ---------------------------------------------------------------
View(Best_SVM_data)

data <- cbind(mydata$status, Best_SVM_data)
colnames(data)[1] <- "status"
# 将数据集分成训练集和测试集（适用于样本量较大的情况）
set.seed(12345)
train_index <- sample(1:nrow(data), nrow(data) * 0.7)# 训练集和测试集7 3分
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# 拟合模型
model_linear <- svm(status~.,data = train_data,kernel="linear",probability = TRUE)

# 查看模型预测准确率
# 模型准确率的计算 formula 为:
#  Accuracy = (TP + TN) / (TP + FP + FN + TN)
# - TP(True Positive):真正例,模型预测为正例,实际也是正例
# - TN(True Negative):真负例,模型预测为负例,实际也是负例
# - FP(False Positive):假正例,模型预测为正例,实际是负例
# - FN(False Negative):假负例,模型预测为负例,实际是正例
mean(train_data[,1] == model_linear$fitted)

# 查看混淆矩阵
table(actual = train_data[,1],model_linear$fitted) 

# 获取模型预测概率
prob_linear <- predict(model_linear, train_data, probability=TRUE)
prob_linear_use <- attr(prob_linear,"probabilities")
linear_ROC <- roc(response = train_data$status,predictor = prob_linear_use[,2])

# plot函数绘制ROC曲线
plot(linear_ROC,
     legacy.axes = TRUE,
     main="ROC curve",
     type = "l",col="red",lty=1,
     print.auc = T,
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best") # 在roc曲线上显示最佳阈值点)




# 拟合四种不同核函数的SVM模型 ---------------------------------------------------------
summary(data)
model_linear <- svm(status~.,data = train_data,kernel="linear",probability = TRUE)
model_polynomial <- svm(status~.,data = train_data,kernel="polynomial",probability = TRUE)
model_radial <- svm(status~.,data = train_data,kernel="radial",probability = TRUE)
model_sigmoid <- svm(status~.,data = train_data,kernel="sigmoid",probability = TRUE)
# 查看模型预测准确率
mean(train_data[,1] == model_linear$fitted)
mean(train_data[,1] == model_polynomial$fitted)
mean(train_data[,1] == model_radial$fitted)
mean(train_data[,1] == model_sigmoid$fitted)
# 查看混淆矩阵
table(actual = train_data[,1],model_linear$fitted) 
table(actual = train_data[,1],model_polynomial$fitted)
table(actual = train_data[,1],model_radial$fitted)
table(actual = train_data[,1],model_sigmoid$fitted)

# 柱状图可视化 ------------------------------------------------------------------
accuracy <- c(mean(train_data[,1] == model_linear$fitted), 
              mean(train_data[,1] == model_polynomial$fitted),
              mean(train_data[,1] == model_radial$fitted),
              mean(train_data[,1] == model_sigmoid$fitted))
mycolors <- c('#54D9ED','#324B4F','#95B0B5','#D1BCFE','#FCFCD4','#F78ABA',
              '#00A9CB','#8EF8B6','#ED546C','#1BBEB9','#DE90DE','#00A44F')
barplot(accuracy, names.arg=c("Linear", "Polynomial", "Radial", "Sigmoid"), 
        main = "Accuracy of model", xlab = "Kernel Type", ylab = "Accuracy",
        col = mycolors)


# 四个模型一起来 -----------------------------------------------------------------
# 获取模型预测概率
prob_linear <- predict(model_linear, train_data, probability=TRUE)  
prob_linear_use <- attr(prob_linear,"probabilities")

prob_polynomial <- predict(model_polynomial, train_data, probability=TRUE)  
prob_polynomial_use <- attr(prob_polynomial,"probabilities")

prob_radial <- predict(model_radial, train_data, probability=TRUE)
prob_radial_use <- attr(prob_radial,"probabilities")  

prob_sigmoid <- predict(model_sigmoid, train_data, probability=TRUE)  
prob_sigmoid_use <- attr(prob_sigmoid,"probabilities")

# 计算ROC曲线  
linear_ROC <- roc(response = train_data$status, predictor = prob_linear_use[,2])  
polynomial_ROC <- roc(response = train_data$status, predictor = prob_polynomial_use[,2])
radial_ROC <- roc(response = train_data$status, predictor = prob_radial_use[,2])   
sigmoid_ROC <- roc(response = train_data$status, predictor = prob_sigmoid_use[,2])

# 提取ROC点  
linear_TPR <- linear_ROC$sensitivities  
linear_FPR <- 1-linear_ROC$specificities

polynomial_TPR <- polynomial_ROC$sensitivities 
polynomial_FPR <- 1-polynomial_ROC$specificities

radial_TPR <- radial_ROC$sensitivities
radial_FPR <- 1-radial_ROC$specificities

sigmoid_TPR <- sigmoid_ROC$sensitivities  
sigmoid_FPR <- 1-sigmoid_ROC$specificities  

# 创建包含FPR和TPR的数据框，用于下一步绘图
ROC_df <- data.frame(FPR = c(linear_FPR, polynomial_FPR, radial_FPR, sigmoid_FPR), 
                     TPR = c(linear_TPR, polynomial_TPR, radial_TPR, sigmoid_TPR),
                     model = c(rep("Linear", length(linear_TPR)), 
                               rep("Polynomial", length(polynomial_TPR)),
                               rep("Radial", length(radial_TPR)),
                               rep("Sigmoid", length(sigmoid_TPR))))

# 绘制ROC曲线 
plot <- ggplot(ROC_df, aes(x = FPR, y = TPR, col = model)) + 
  geom_line(size = 0.8) +
  geom_line(color = "#BC1328", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Linear")) +     # Linear模型,红色线  
  geom_line(color = "#0067B5", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Polynomial")) + # Polynomial模型,蓝色线
  geom_line(color = "#09891D", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Radial")) +   # Radial模型,绿色线   
  geom_line(color = "orange", linetype = 1, size = 0.8, data = subset(ROC_df, model == "Sigmoid")) + # Sigmoid模型,橙色线
  # 添加对角线
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +  
  # 设置标题和轴
  ggtitle("ROC Curves") +    
  labs(x = "1-specificity", y = "Sensitivity") +
  # 修改轴文本和标题
  theme(axis.text = element_text(face = "bold", size = 14, color = "black"),    
        axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))+
  # 设置图例位置
  theme(legend.position="bottom", legend.box = "horizontal")


#设置不同主题风格
plot + theme_grey() #灰色主题
plot + theme_bw()  #黑白主题
plot + theme_classic() #经典主题