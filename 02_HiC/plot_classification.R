library(ggplot2)
library(scales)

projdir = "C:/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C/02_HIC"

# young
# total inter-chromosomal: 43486741 (23.3894%)
# total intra-chromosomal: 142438325 (76.6106%)
# short-range: 26122731 (18.3397% / 14.0501%)
# long-range: 116315594 (81.6603% / 62.5605%)

# aged
# total inter-chromosomal: 34516118 (22.0723%)
# total intra-chromosomal: 121861375 (77.9277%)
# short-range: 22102408 (18.1373% / 14.134%)
# long-range: 99758967 (81.8627% / 63.7937%)

df = data.frame(matrix(c(43486741,142438325,26122731,116315594,
                34516118,121861375,22102408,99758967),
              nrow = 2, byrow = T))
colnames(df) = c("Total inter-chromosomal", "Total intra-chromosomal", "Short-cis (<20kb)", "Long-cis (>20kb)")
df$Age = c("Y","A")

plt.df = df %>%
  pivot_longer(cols = 1:4, 
               names_to="Interactions") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Age = factor(Age, levels=c("Y","A"))) %>%
  group_by(Age) %>%
  mutate(Fraction = value/(value[Interactions=="Total inter-chromosomal"] + value[Interactions=="Total intra-chromosomal"]))

plt.df %>%
  filter(Interactions!="Total intra-chromosomal") %>%
  ggplot(aes(x="", y=value, fill=Interactions)) +
  geom_col(width=1, position = position_fill()) +
  geom_text(aes(label = sprintf("%s\n(%s)", scientific(value), percent(Fraction,accuracy = 0.1))), size=3, position = position_fill(vjust = 0.5)) +
  facet_wrap(~Age) +
  coord_polar("y", start=0) +
  scale_y_continuous(labels = scales::percent) +
  xlab("") + ylab("") +
  theme_pubclean() + 
  scale_fill_brewer(palette="Blues") +
  theme(axis.text = element_text(size=8, color="black"),
        legend.text = element_text(size=10),
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow=3))
ggsave(file.path(projdir, "interaction_classification.png"), dpi=300, width=5, height=5)

plt.df %>%
  group_by(Age) %>%
  summarise(short.prop = value[name=="short-range"]/value[name=="total intra-chromosomal"],
            long.prop = value[name=="long-range"]/value[name=="total intra-chromosomal"])
