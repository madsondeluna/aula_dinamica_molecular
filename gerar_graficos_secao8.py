#!/usr/bin/env python3
"""
Script para gerar gráficos específicos para a seção 8 do README
"""

import matplotlib.pyplot as plt
import numpy as np
import os

# Configurações
plt.style.use('default')
np.random.seed(42)
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12

def create_rmsd_plot():
    """Cria gráfico de exemplo para RMSD"""
    time = np.linspace(0, 10, 1000)
    rmsd = 0.15 + 0.05 * np.random.random(1000) + 0.02 * np.sin(time * 2)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, rmsd, color='#2E8B57', linewidth=2, alpha=0.8)
    plt.fill_between(time, rmsd, alpha=0.3, color='#2E8B57')
    
    # Linha da média
    mean_rmsd = np.mean(rmsd)
    plt.axhline(y=mean_rmsd, color='red', linestyle='--', linewidth=2, 
                label=f'Média: {mean_rmsd:.3f} nm')
    
    plt.xlabel('Tempo (ns)', fontsize=14, fontweight='bold')
    plt.ylabel('RMSD (nm)', fontsize=14, fontweight='bold')
    plt.title('RMSD do Backbone da Lisozima', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('imagens/rmsd_exemplo.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_rmsf_plot():
    """Cria gráfico de exemplo para RMSF"""
    residues = np.arange(1, 130)
    rmsf = 0.1 + 0.15 * np.random.random(129)
    # Simula loops flexíveis
    rmsf[10:15] *= 2.5
    rmsf[45:50] *= 2.0
    rmsf[80:85] *= 2.3
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(residues, rmsf, width=1.0, color='#4169E1', alpha=0.7, edgecolor='black', linewidth=0.3)
    
    # Destaca resíduos flexíveis
    for i, rmsf_val in enumerate(rmsf):
        if rmsf_val > 0.25:
            bars[i].set_color('#FF6347')
    
    plt.axhline(y=np.mean(rmsf), color='green', linestyle='--', linewidth=2,
                label=f'Média: {np.mean(rmsf):.3f} nm')
    plt.axhline(y=0.25, color='red', linestyle=':', linewidth=2,
                label='Limite alta flexibilidade (0.25 nm)')
    
    plt.xlabel('Número do Resíduo', fontsize=14, fontweight='bold')
    plt.ylabel('RMSF (nm)', fontsize=14, fontweight='bold')
    plt.title('Flutuação por Resíduo (C-alpha)', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='y')
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('imagens/rmsf_exemplo.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_gyration_plot():
    """Cria gráfico de exemplo para raio de giração"""
    time = np.linspace(0, 10, 1000)
    gyration = 1.45 + 0.02 * np.random.random(1000) + 0.005 * np.sin(time * 3)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, gyration, color='#FF8C00', linewidth=2)
    plt.fill_between(time, gyration, alpha=0.3, color='#FF8C00')
    
    mean_gyration = np.mean(gyration)
    std_gyration = np.std(gyration)
    
    plt.axhline(y=mean_gyration, color='blue', linestyle='--', linewidth=2,
                label=f'Média: {mean_gyration:.3f} nm')
    
    plt.fill_between(time, 
                    mean_gyration - std_gyration, 
                    mean_gyration + std_gyration, 
                    alpha=0.2, color='blue', label=f'±1σ ({std_gyration:.3f} nm)')
    
    plt.xlabel('Tempo (ns)', fontsize=14, fontweight='bold')
    plt.ylabel('Raio de Giração (nm)', fontsize=14, fontweight='bold')
    plt.title('Raio de Giração da Proteína', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('imagens/giracao_exemplo.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_hbond_plot():
    """Cria gráfico de exemplo para ligações de hidrogênio"""
    time = np.linspace(0, 10, 1000)
    hbonds = 115 + 10 * np.random.random(1000) + 3 * np.sin(time * 4)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, hbonds, color='#9932CC', linewidth=2, alpha=0.8)
    plt.fill_between(time, hbonds, alpha=0.3, color='#9932CC')
    
    mean_hbond = np.mean(hbonds)
    plt.axhline(y=mean_hbond, color='red', linestyle='--', linewidth=2,
                label=f'Média: {mean_hbond:.1f} ligações')
    
    # Média móvel
    window_size = 50
    moving_avg = np.convolve(hbonds, np.ones(window_size)/window_size, mode='valid')
    time_ma = time[window_size-1:]
    plt.plot(time_ma, moving_avg, color='black', linewidth=2, 
             label=f'Média móvel (janela: {window_size})')
    
    plt.xlabel('Tempo (ns)', fontsize=14, fontweight='bold')
    plt.ylabel('Número de Ligações H', fontsize=14, fontweight='bold')
    plt.title('Ligações de Hidrogênio Intramoleculares', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('imagens/hbond_exemplo.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    print("Gerando gráficos para a seção 8 do README...")
    
    # Cria pasta de imagens se não existir
    if not os.path.exists('imagens'):
        os.makedirs('imagens')
    
    create_rmsd_plot()
    print("✓ Gráfico RMSD criado: imagens/rmsd_exemplo.png")
    
    create_rmsf_plot()
    print("✓ Gráfico RMSF criado: imagens/rmsf_exemplo.png")
    
    create_gyration_plot()
    print("✓ Gráfico Raio de Giração criado: imagens/giracao_exemplo.png")
    
    create_hbond_plot()
    print("✓ Gráfico Ligações H criado: imagens/hbond_exemplo.png")
    
    print("\nTodos os gráficos foram gerados com sucesso!")

if __name__ == "__main__":
    main()
