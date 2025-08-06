#!/usr/bin/env python3
"""
Script para gerar gráficos de exemplo para análise de dinâmica molecular
"""

import matplotlib.pyplot as plt
import numpy as np
import os

# Configurações gerais
plt.style.use('seaborn-v0_8')
np.random.seed(42)

def generate_rmsd_plot():
    """Gera gráfico de exemplo para RMSD"""
    time = np.linspace(0, 10, 1000)
    rmsd = 0.15 + 0.05 * np.random.random(1000) + 0.02 * np.sin(time * 2)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, rmsd, color='#2E8B57', linewidth=1.5, alpha=0.8)
    plt.fill_between(time, rmsd, alpha=0.3, color='#2E8B57')
    
    mean_rmsd = np.mean(rmsd)
    plt.axhline(y=mean_rmsd, color='red', linestyle='--', linewidth=2, 
                label=f'Média: {mean_rmsd:.3f} nm')
    
    plt.xlabel('Tempo (ns)', fontsize=14)
    plt.ylabel('RMSD (nm)', fontsize=14)
    plt.title('RMSD do Backbone da Lisozima', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    # Adiciona estatísticas
    textstr = f'Desvio padrão: {np.std(rmsd):.3f} nm\\nValor final: {rmsd[-1]:.3f} nm'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('exemplo_rmsd.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_rmsf_plot():
    """Gera gráfico de exemplo para RMSF"""
    residues = np.arange(1, 130)
    rmsf = 0.1 + 0.15 * np.random.random(129)
    # Simula loops mais flexíveis
    rmsf[10:15] *= 2
    rmsf[45:50] *= 1.8
    rmsf[80:85] *= 2.2
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(residues, rmsf, width=0.8, color='#4169E1', alpha=0.7, 
                   edgecolor='black', linewidth=0.5)
    
    # Destaca resíduos flexíveis
    for i, (res, rmsf_val) in enumerate(zip(residues, rmsf)):
        if rmsf_val > 0.25:
            bars[i].set_color('#FF6347')
    
    plt.axhline(y=np.mean(rmsf), color='green', linestyle='--', linewidth=2,
                label=f'Média: {np.mean(rmsf):.3f} nm')
    plt.axhline(y=0.25, color='red', linestyle=':', linewidth=2,
                label='Limite alta flexibilidade (0.25 nm)')
    
    plt.xlabel('Número do Resíduo', fontsize=14)
    plt.ylabel('RMSF (nm)', fontsize=14)
    plt.title('Flutuação por Resíduo (C-alpha)', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='y')
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    plt.savefig('exemplo_rmsf.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_gyration_plot():
    """Gera gráfico de exemplo para raio de giração"""
    time = np.linspace(0, 10, 1000)
    gyration = 1.45 + 0.02 * np.random.random(1000) + 0.005 * np.sin(time * 3)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, gyration, color='#FF8C00', linewidth=1.5)
    plt.fill_between(time, gyration, alpha=0.3, color='#FF8C00')
    
    mean_gyration = np.mean(gyration)
    std_gyration = np.std(gyration)
    
    plt.axhline(y=mean_gyration, color='blue', linestyle='--', linewidth=2,
                label=f'Média: {mean_gyration:.3f} nm')
    
    plt.fill_between(time, 
                    mean_gyration - std_gyration, 
                    mean_gyration + std_gyration, 
                    alpha=0.2, color='blue', label=f'±1σ ({std_gyration:.3f} nm)')
    
    plt.xlabel('Tempo (ns)', fontsize=14)
    plt.ylabel('Raio de Giração (nm)', fontsize=14)
    plt.title('Raio de Giração da Proteína', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    variation_percent = (std_gyration / mean_gyration) * 100
    textstr = f'Variação: {variation_percent:.2f}%\\nStabilidade: {"Alta" if variation_percent < 2 else "Moderada" if variation_percent < 5 else "Baixa"}'
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('exemplo_raio_giracao.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_hbond_plot():
    """Gera gráfico de exemplo para ligações de hidrogênio"""
    time = np.linspace(0, 10, 1000)
    hbonds = 115 + 10 * np.random.random(1000) + 3 * np.sin(time * 4)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, hbonds, color='#9932CC', linewidth=1.5, alpha=0.8)
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
    
    plt.xlabel('Tempo (ns)', fontsize=14)
    plt.ylabel('Número de Ligações H', fontsize=14)
    plt.title('Ligações de Hidrogênio Intramoleculares', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    coef_var = (np.std(hbonds)/mean_hbond)*100
    textstr = f'Coef. Variação: {coef_var:.1f}%\\nEstabilidade: {"Alta" if coef_var < 10 else "Moderada" if coef_var < 15 else "Baixa"}'
    props = dict(boxstyle='round', facecolor='plum', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('exemplo_ligacoes_h.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_combined_plot():
    """Gera gráfico combinado com todas as análises"""
    # Dados simulados
    time = np.linspace(0, 10, 1000)
    rmsd = 0.15 + 0.05 * np.random.random(1000) + 0.02 * np.sin(time * 2)
    gyration = 1.45 + 0.02 * np.random.random(1000) + 0.005 * np.sin(time * 3)
    hbonds = 115 + 10 * np.random.random(1000) + 3 * np.sin(time * 4)
    
    residues = np.arange(1, 130)
    rmsf = 0.1 + 0.15 * np.random.random(129)
    rmsf[10:15] *= 2
    rmsf[45:50] *= 1.8
    rmsf[80:85] *= 2.2
    
    # Cria figura combinada
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Análise Completa de Dinâmica Molecular - Lisozima', 
                 fontsize=18, fontweight='bold')
    
    # RMSD
    axes[0,0].plot(time, rmsd, color='#2E8B57', linewidth=1.5)
    axes[0,0].axhline(y=np.mean(rmsd), color='red', linestyle='--', alpha=0.8)
    axes[0,0].set_xlabel('Tempo (ns)')
    axes[0,0].set_ylabel('RMSD (nm)')
    axes[0,0].set_title('RMSD do Backbone')
    axes[0,0].grid(True, alpha=0.3)
    
    # RMSF
    bars = axes[0,1].bar(residues, rmsf, width=0.8, color='#4169E1', alpha=0.7)
    for i, rmsf_val in enumerate(rmsf):
        if rmsf_val > 0.25:
            bars[i].set_color('#FF6347')
    axes[0,1].axhline(y=np.mean(rmsf), color='green', linestyle='--', alpha=0.8)
    axes[0,1].set_xlabel('Número do Resíduo')
    axes[0,1].set_ylabel('RMSF (nm)')
    axes[0,1].set_title('Flutuação por Resíduo')
    axes[0,1].grid(True, alpha=0.3)
    
    # Raio de Giração
    axes[1,0].plot(time, gyration, color='#FF8C00', linewidth=1.5)
    axes[1,0].axhline(y=np.mean(gyration), color='blue', linestyle='--', alpha=0.8)
    axes[1,0].set_xlabel('Tempo (ns)')
    axes[1,0].set_ylabel('Raio de Giração (nm)')
    axes[1,0].set_title('Raio de Giração')
    axes[1,0].grid(True, alpha=0.3)
    
    # Ligações de Hidrogênio
    axes[1,1].plot(time, hbonds, color='#9932CC', linewidth=1.5)
    axes[1,1].axhline(y=np.mean(hbonds), color='red', linestyle='--', alpha=0.8)
    axes[1,1].set_xlabel('Tempo (ns)')
    axes[1,1].set_ylabel('Número de Ligações H')
    axes[1,1].set_title('Ligações de Hidrogênio')
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('exemplo_analise_completa.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Função principal que gera todos os gráficos"""
    print("Gerando gráficos de exemplo para análise de dinâmica molecular...")
    
    # Cria pasta para as imagens se não existir
    if not os.path.exists('imagens_exemplos'):
        os.makedirs('imagens_exemplos')
    
    os.chdir('imagens_exemplos')
    
    # Gera todos os gráficos
    generate_rmsd_plot()
    print("✓ Gráfico RMSD gerado")
    
    generate_rmsf_plot()
    print("✓ Gráfico RMSF gerado")
    
    generate_gyration_plot()
    print("✓ Gráfico Raio de Giração gerado")
    
    generate_hbond_plot()
    print("✓ Gráfico Ligações de Hidrogênio gerado")
    
    generate_combined_plot()
    print("✓ Gráfico Combinado gerado")
    
    print(f"\\nTodos os gráficos foram salvos na pasta 'imagens_exemplos'")
    print("Arquivos gerados:")
    print("- exemplo_rmsd.png")
    print("- exemplo_rmsf.png")
    print("- exemplo_raio_giracao.png")
    print("- exemplo_ligacoes_h.png")
    print("- exemplo_analise_completa.png")

if __name__ == "__main__":
    main()
