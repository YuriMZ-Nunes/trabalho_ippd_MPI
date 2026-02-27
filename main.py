import pandas as pd
import matplotlib.pyplot as plt
import glob
import imageio.v2 as imageio
import os
import numpy as np

print("Lendo as planilhas da simulação (Agentes e Território)...")
df_agentes = pd.concat([pd.read_csv(f) for f in glob.glob('agentes_rank_*.csv')], ignore_index=True)
df_grid = pd.concat([pd.read_csv(f) for f in glob.glob('grid_rank_*.csv')], ignore_index=True)

ciclos = sorted(df_agentes['ciclo'].unique())
frames = []

# Descobre o tamanho do mapa pelas coordenadas máximas salvas
W_max = df_grid['x'].max() + 1
H_max = df_grid['y'].max() + 1

print("Pintando os quadros (frames) da animação...")
for c in ciclos:
    agentes_c = df_agentes[df_agentes['ciclo'] == c]
    grid_c = df_grid[df_grid['ciclo'] == c]
    
    # Monta a matriz do mapa de calor usando NumPy para ser rápido
    matriz_recursos = np.zeros((W_max, H_max))
    matriz_recursos[grid_c['x'].values, grid_c['y'].values] = grid_c['recurso'].values
    
    plt.figure(figsize=(8, 8))
    plt.title(f"Mobilidade Sazonal - Ciclo {c:03d}")
    
    # 1. Pinta o Mapa de Fundo (Transposto para alinhar X e Y corretamente)
    # cmap='YlGn' vai do Amarelo (seco/vazio) para o Verde Escuro (cheio de recursos)
    plt.imshow(matriz_recursos.T, cmap='YlGn', origin='lower', vmin=0, vmax=150)
    
    # 2. Plota os Agentes por cima (Pontos vermelhos)
    plt.scatter(agentes_c['x'], agentes_c['y'], c='red', s=10, edgecolor='black', linewidth=0.5)
    
    plt.xlim(0, W_max - 1)
    plt.ylim(0, H_max - 1)
    
    nome_imagem = f"frame_{c:03d}.png"
    plt.savefig(nome_imagem, dpi=100)
    plt.close()
    
    frames.append(imageio.imread(nome_imagem))

print("Costurando os quadros em um GIF fluido...")
# fps=15 deixa a animação bem fluida (15 quadros por segundo)
imageio.mimsave('animacao_hibrida.gif', frames, fps=15)

# Limpeza das imagens PNG soltas
for c in ciclos:
    os.remove(f"frame_{c:03d}.png")

print("Pronto! Abra o 'animacao_hibrida.gif' para ver o resultado.")