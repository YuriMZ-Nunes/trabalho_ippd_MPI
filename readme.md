# Simulação de Mobilidade Sazonal Híbrida (MPI + OpenMP)

Este projeto implementa uma simulação complexa de agentes em um ambiente dinâmico com recursos sazonais. A aplicação utiliza uma arquitetura híbrida de computação paralela, combinando **MPI** para a decomposição de domínio (distribuição entre diferentes processos/nós) e **OpenMP** para o paralelismo de grão fino dentro de cada processo.



## 🚀 Funcionalidades

* **Arquitetura Híbrida**: Divisão do território em subgrids (MPI) com processamento acelerado de agentes via threads (OpenMP).
* **Troca de Halos (Ghost Cells)**: Sincronização de bordas entre processos para permitir que agentes tomem decisões baseadas em células vizinhas de outros ranks.
* **Dinâmica Sazonal**: O ambiente alterna entre estações (Seca/Cheia), afetando a taxa de regeneração dos recursos no grid.
* **Comunicação Não-Bloqueante**: Uso de `MPI_Isend` e `MPI_Irecv` para otimizar a migração de agentes e a atualização de halos sem travar o processamento.
* **Visualização**: Script Python integrado para gerar uma animação `.gif` a partir dos logs de simulação.

## 🛠️ Tecnologias Utilizadas

* **C++**: Núcleo da simulação.
* **MPI (Message Passing Interface)**: Coordenação e distribuição de carga entre processos.
* **OpenMP**: Paralelismo multi-core para movimentação de agentes e atualização de recursos.
* **Python (Pandas/Matplotlib/ImageIO)**: Pós-processamento e geração de visualizações.

## 🏗️ Estrutura do Projeto

* `main.cpp`: Código-fonte principal com a lógica de simulação, tipos derivados MPI e kernels OpenMP.
* `main.py`: Script para leitura dos arquivos `.csv` e criação da animação.
* `makefile`: Automação da compilação, execução e limpeza do ambiente.

## ⚙️ Como Executar

### Pré-requisitos
Certifique-se de ter instalado:
* Compilador de MPI (ex: `mpich` ou `openmpi`).
* Bibliotecas OpenMP.
* Python 3 com `pandas`, `matplotlib` e `imageio`.

### Passo a passo
1.  **Compilar e Rodar**:
    O comando abaixo compila o código, executa a simulação com 4 processos e gera o GIF automaticamente:
    ```bash
    make run
    ```

2.  **Limpar Arquivos**:
    Para remover executáveis e logs gerados:
    ```bash
    make clean
    ```

## 📊 Detalhes Técnicos

### Decomposição de Domínio
O mundo de $200 \times 200$ é dividido verticalmente entre os processos MPI. Cada processo é responsável por uma fatia local de largura `local_W = W / size`.

### Tipos Derivados MPI
Para otimizar a comunicação, foram criados tipos estruturados (`MPI_Type_create_struct`) para as structs `Agent` e `Cell`, permitindo o envio direto de buffers de memória sem serialização manual.

### Equilíbrio de Carga
A função `executar_carga(r)` simula um processamento computacional proporcional à quantidade de recursos na célula, desafiando o escalonamento das threads OpenMP.

---
