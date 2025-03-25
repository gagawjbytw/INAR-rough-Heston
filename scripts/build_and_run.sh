#!/bin/bash

# 切换到项目根目录
cd "$(dirname "$0")/.."

echo "Compiling monte_carlo_inar_all_options.cpp..."

# 编译源文件
g++ -std=c++17 src/monte_carlo_inar_all_options.cpp -o monte_carlo_inar_all_options

# 检查编译是否成功
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    echo "Please make sure you have g++ installed"
    read -p "Press Enter to continue..."
    exit 1
fi

echo "Compilation successful!"
echo
echo "Running the program..."
echo

# 运行编译后的程序
./monte_carlo_inar_all_options

echo
read -p "Press Enter to continue..." 
