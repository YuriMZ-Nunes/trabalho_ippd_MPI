# Define o compilador wrapper do MPI para C++
CXX = mpicxx

# Flags de compilação: 
# -O3: Otimização máxima de desempenho
# -Wall: Mostra todos os avisos (warnings)
# -fopenmp: Habilita as diretivas do OpenMP essenciais para o paralelismo intra-nó
CXXFLAGS = -O3 -Wall -fopenmp

# Nome do arquivo executável de saída
TARGET = simulacao

# Código fonte principal
SRCS = main.cpp

# Regra principal que é executada ao digitar apenas 'make'
all: $(TARGET)

# Regra para gerar o executável
$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

# Regra para limpar os binários, dados da simulação e imagens (digite 'make clean')
clean:
	rm -f $(TARGET) *.csv *.png *.gif

# Regra completa: compila, limpa dados antigos, roda a simulação híbrida e gera o GIF
run: $(TARGET)
	@echo "--- 1. Limpando dados antigos ---"
	rm -f *.csv
	@echo "--- 2. Iniciando Simulação MPI + OpenMP ---"
	mpirun -np 4 ./$(TARGET)
	@echo "--- 3. Iniciando Renderização em Python ---"
	python3 main.py