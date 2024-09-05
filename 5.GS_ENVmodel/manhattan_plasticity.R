library(getopt,quietly = T)
command=matrix(c( 
  'help', 'h', 0,'loical', 'help information',
  'map', 'M', 1, 'character', 'Physical coordinate of Markers , Contains three columns, namely SNP CHR POS.',
  'snplist', 'S', 1, 'character', 'SNP Contribution Table, Contains two columns, IID featureGain_sum.',
  'mingain','G',1,'num','The coefficient that controls the threshold is obtained by multiplying the maximum contribution value of a single SNP by this coefficient'
),
byrow=T,ncol=5)
args=getopt(command)

if (!is.null(args$help) || is.null(args$map) || is.null(args$snplist) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

## 设置默认值
if ( is.null(args$mingain) ) {
 args$mingain = 0.05
}

library(readxl,quietly = T)
library(dplyr,quietly = T)
library(qqman,quietly = T)
library(ggplot2)
library(tidyverse)

map_file <- args$map
snplist_file <- args$snplist

map <- read_excel(map_file, col_names = T)
snp <- read.table(snplist_file, header=F, sep=",")
colnames(snp) <- c("SNP", "FeatureGainSum", "group")
group <- unique(snp$group)

data=NULL

for ( i in group){
  snp.group <- snp[which(snp$group==i),]
  # 合并 map 和 snp 数据框
  data.group <- merge(map, snp.group, by = "SNP", all.x = TRUE)

  # 填充特征值为 0，把不完整的分组补全
  data.group$FeatureGainSum <- ifelse(is.na(data.group$FeatureGainSum), 0, data.group$FeatureGainSum)
  data.group$group <- ifelse(is.na(data.group$group), i, data.group$group)

  data <- rbind(data, data.group)

}



# 1)计算chr长度
chr_len <- data %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(POS))
# 2） 计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)
#3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(data, ., by="CHR") %>%
  arrange(CHR, POS) %>%
  mutate( POScum = POS + total)

#计算每个snp累计贡献和单个snp最大贡献的比例，选择0.05以上的
max_feature_gain <- Snp_pos %>%
  group_by(group) %>%
  summarise(max_feature_gain = max(FeatureGainSum))

# 合并最大 FeatureGainSum 到原始数据框
Snp_pos <- Snp_pos %>%
  left_join(max_feature_gain, by = "group") %>%
  # 计算 ratio
  mutate(ratio = FeatureGainSum / max_feature_gain)

# 选择满足条件的 SNP
select_snp <- Snp_pos %>%
  filter(ratio > args$mingain)

# 将结果转换为数据框
select_snp <- as.data.frame(select_snp)

#去除重复SNP
select_snp_filtered <- select_snp %>%  
  group_by(SNP) %>%  
  filter(ratio == max(ratio)) %>%  
  ungroup()  

# 将结果转换为数据框
select_snp_filtered <- as.data.frame(select_snp_filtered)

colnames(select_snp)[1] <- "IID"
colnames(select_snp_filtered)[1] <- "IID"

X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(POScum) +min(POScum) ) / 2 )

color <- c("#367db8","#fb9a9b","#e4191c","#fcb359","#87d5c6","#fecbe3","#be7fbc","#ab5326","#49b243","#9fb2da","#89a6af")

# 定义图例标签
legend_labels <- unique(Snp_pos$group)

p <- ggplot(Snp_pos, aes(x = POScum, y = ratio)) +
  geom_point(aes(color = as.factor(group)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = color, labels = legend_labels) +  # 添加图例标签
  scale_x_continuous(breaks = X_axis$center, labels = X_axis$CHR) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Ratio to the max feature of SNP", color = "Parameter type") + 
#  geom_hline(yintercept = args$mingain, color = 'black', size = 0.6, linetype = 'twodash') +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 2, r = 1, b = 1, l = 1, unit = "cm") # 设置上边界距离为2cm
  ) +
  coord_cartesian(ylim = c(0, 1.1)) 
# + geom_text(data = select_snp_filtered, aes(x = POScum, y = ratio, label = IID), vjust = -0.5, size = 3)

# 调整图例位置
# p + theme(legend.position = "top")

ggsave("Plasticity_select_snp_feature_manhattan.png", plot = p, width = 9, height = 6, dpi = 600)

select_snp <- as.data.frame(select_snp[,1])
colnames(select_snp)[1] <- "IID"
write.table(select_snp, file = "Plasticity_select_snp.list", col.names = F, row.names = F, quote = F )