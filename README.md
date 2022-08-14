# IPCC2022-Pivot

## 编译 
使用makefile，自动创建build目录并将二进制可执行文件pivot编译并输出至./build目录


用-i选项忽略掉build目录已存在的情况
```
make -i
```
## 运行

- submit.sh中指定了运行参数与执行命令


- run.sh通过sbatch提交submit.sh

```
./run.sh $INPUT

./run.sh uniformvector-2dim-5h.txt

./run.sh uniformvector-4dim-1h.txt

```
