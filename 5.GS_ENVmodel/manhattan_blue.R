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
snp <- read.table(snplist_file, header=T, sep=",")
colnames(snp) <- c("SNP", "FeatureGainSum")

# 合并 map 和 snp 数据框
data <- merge(map, snp, by = "SNP", all.x = TRUE)

# 填充缺失值为 0
data$FeatureGainSum <- ifelse(is.na(data$FeatureGainSum), 0, data$FeatureGainSum)

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
Snp_pos$ratio <- Snp_pos$FeatureGainSum/max(Snp_pos$FeatureGainSum)
select_snp <- as.data.frame(Snp_pos[which(Snp_pos$ratio>args$mingain),])
colnames(select_snp)[1] <- "IID"

X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(POScum) +min(POScum) ) / 2 )

p <- ggplot(Snp_pos, aes(x = POScum, y = ratio)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep("#367db8", length(unique(Snp_pos$CHR)))) +
  scale_x_continuous(breaks = X_axis$center, labels = X_axis$CHR) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Ratio to the max feature of SNP") + 
#  geom_hline(yintercept = args$mingain, color = 'black', size = 0.6, linetype = 'twodash') +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 2, r = 1, b = 1, l = 1, unit = "cm") # 设置上边界距离为2cm
  ) +
  coord_cartesian(ylim = c(0, 1.1)) 
#  +geom_text(data = select_snp, aes(x = POScum, y = ratio, label = IID), vjust = -0.5, size = 3) # 添加文本标签

ggsave("BLUP_select_snp_feature_manhattan.png", plot = p, width = 9, height = 6, dpi = 600)

select_snp <- as.data.frame(select_snp[,1])
colnames(select_snp)[1] <- "IID"
write.table(select_snp, file = "BLUP_select_snp.list", col.names = T, row.names = F, quote = F )