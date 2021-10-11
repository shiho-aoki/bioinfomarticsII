#
# bioinfomatices II 
# Day2 遺伝子発現解析
#
# 2021/10/11 Shiho AOKI

# install packages.
# install.packages("tidyverse");
# install.packages("ggfortify");

library(tidyverse);
library(ggfortify);
rm(list=ls());

# genoID: genome name
# cntx: expression level of control group
# trtx: expression level of treated group
GENO_EXP_DATA <- read.delim("/Users/shihoaoki/main/connectGoogle/gene_expression.txt", header = TRUE);

# observation of data
summary(GENO_EXP_DATA);
# check boxplot
boxplot(GENO_EXP_DATA[, 2:11]);
boxplot(log(GENO_EXP_DATA[, 2:11]));
# check boxplot drawn by ggplot (exclude geneID)
GENO_EXP_DATA_LONG <- gather(GENO_EXP_DATA[, 2:11], "Condition", "Expression");
ggplot(GENO_EXP_DATA_LONG, aes(x=Condition, y=Expression)) + geom_boxplot();
ggplot(GENO_EXP_DATA_LONG, aes(x=Condition, y=log(Expression))) + geom_boxplot();

# 発現量の対数値をt検定する
# →i番目とj番目の遺伝子の発現量の平均値が有意に違うかどうかを検定する
LogDATA <- log(GENO_EXP_DATA[, 2:11]);
i <- 1
UP_LIST <- c();
N <- length(LogDATA[,1]);
p <- rep(0, N); 
for (i in 1:N) {
  p[i] <- t.test(LogDATA[i, 1:5],
                 LogDATA[i, 6:10])$p.value;
  # p[i]が0.05以下のものは発現変動遺伝子とみなす
  if(p[i] < 0.05 && mean(as.matrix(LogDATA[i,1:5])) < mean(as.matrix(LogDATA[i,6:10])) ){
    UP_LIST <- c(UP_LIST, GENO_EXP_DATA[i,"geneID"]);
  };
}
hist(p);

# 発現変動遺伝子を調べる
GO_TERM <- read.delim("/Users/shihoaoki/main/connectGoogle/GOlist13359.tsv", 
                      header = FALSE, stringsAsFactors = FALSE);
TERMS <- gather(GO_TERM[2:12], col, GOterm, 1:11);
TERMS <- sort(TERMS[, 2]);
TERMS <- terms[TERMS != ""];
M <- 1
for (i in 2:length(TERMS)) {
  if (TERMS[i] != TERMS[i-1]) {
    M <- M+1;
  }
}
L <- 1;
GO_TERMS <- rep("", M);
GO_TERMS[1] <- TERMS[1];
for (i in 2:length(TERMS)) {
  if (TERMS[1] != GO_TERMS[L]){
    L <- L + 1; GO_TERMS[L] <- TERMS[i];
  }
  if (L >= M){ break };
};

dim(GO_TERM);
# 発現量データ全体で、各 GO term が何回ずつ出てくるか（それが付いている遺伝子の数）を数える
GO_UNITED <- unite(GO_TERM[, 2:12], "terms");
NUM_GO_TERMS <- rep(0, M);
for(i in 1:length(GO_TERMS)){
  SUBSTR_GO <- substr(GO_TERMS[i], 1, 10);
  NUM_GO_TERMS[i] <- length(grep(SUBSTR_GO, as.matrix(GO_UNITED)));
}
# p値の小遺伝子のアノテーションに各GOTermが何回ずつ出てくるかを数える
P_VAL <- p[1:length(GO_TERM[,1])];
P_VAL_UNITED <- GO_UNITED[P_VAL <= 0.05];
NUM_TERM_SPECIAL_GENE <- rep(0, length(GO_TERMS));
for (i in 1: length(GO_TERMS)) {
  SUBSTR_GO <- substr(GO_TERMS[i], 1, 10);
  if (length(grep(SUBSTR_GO, as.matrix(P_VAL_UNITED)))<1){ next };
  NUM_TERM_SPECIAL_GENE[i] <- length(grep(SUBSTR_GO, as.matrix(P_VAL_UNITED)));
}

for (i in 1:length(GO_TERMS)) {
  x <- NUM_TERM_SPECIAL_GENE[i];
  y <- length(P_VAL[P_VAL <= 0.05]) -x;
  z <- GO_TERMS[i] - x;
  a <- M - z - y - x;
  P_VAL_OF_GO_TERM <- fisher.test(
    matrix(y(a, z, y, x), ncol=2))$p.value;
  if (P_VAL_OF_GO_TERM < 0.05) print(M[i]);
}
