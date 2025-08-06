# Guia de Análise e Visualização de Resultados

Este guia complementa o tutorial principal com informações detalhadas sobre como analisar e visualizar os resultados da simulação de dinâmica molecular da lisozima.

## Ferramentas de Análise

### 1. Scripts Python

#### Script para Gerar Gráficos de Exemplo
```bash
python gerar_graficos_exemplo.py
```

Este script cria:
- `exemplo_rmsd.png` - Gráfico RMSD com estatísticas
- `exemplo_rmsf.png` - Gráfico RMSF destacando resíduos flexíveis
- `exemplo_raio_giracao.png` - Gráfico do raio de giração com análise de estabilidade
- `exemplo_ligacoes_h.png` - Gráfico de ligações de hidrogênio com média móvel
- `exemplo_analise_completa.png` - Painel combinado de todas as análises

#### Notebook Jupyter Interativo
```bash
jupyter notebook analise_graficos_md.ipynb
```

### 2. Dependências Necessárias

Para executar as análises, instale as seguintes bibliotecas Python:

```bash
pip install matplotlib numpy pandas scipy seaborn jupyter
```

## Interpretação dos Resultados

### RMSD (Root-Mean-Square Deviation)

**O que mede:** Desvio médio quadrático da estrutura em relação à configuração inicial.

**Valores esperados para lisozima:**
- Estrutura estável: 0.1-0.3 nm
- Estrutura moderadamente flexível: 0.3-0.5 nm
- Possível desnaturação: >0.5 nm

**Como interpretar:**
- **Plateau inicial:** Indica equilibração do sistema
- **Tendência crescente:** Possível instabilidade ou mudança conformacional
- **Oscilações:** Movimentos térmicos normais
- **Valor final estável:** Boa convergência da simulação

### RMSF (Root-Mean-Square Fluctuation)

**O que mede:** Flexibilidade de cada resíduo da proteína.

**Valores típicos:**
- Estruturas secundárias (α-hélices, folhas-β): <0.2 nm
- Loops e turns: 0.2-0.4 nm
- Terminais N e C: >0.3 nm

**Regiões importantes na lisozima:**
- Resíduos 10-15: Loop flexível próximo ao sítio ativo
- Resíduos 45-50: Loop de superfície
- Resíduos 80-85: Região de ligação ao substrato
- Terminais: Naturalmente flexíveis

### Raio de Giração

**O que mede:** Compactação global da proteína.

**Valores para lisozima:**
- Estrutura compacta: ~1.4-1.6 nm
- Variação normal: <5% do valor médio
- Instabilidade: >10% de variação

**Indicadores de estabilidade:**
- **Baixa variação (<2%):** Estrutura muito estável
- **Variação moderada (2-5%):** Estrutura estável com flexibilidade normal
- **Alta variação (>5%):** Possível instabilidade ou transição conformacional

### Ligações de Hidrogênio

**O que mede:** Número de ligações de hidrogênio intramoleculares.

**Valores para lisozima:**
- Estrutura nativa: ~100-130 ligações
- Flutuação normal: ±10-15 ligações
- Perda significativa (>20%): Indica desnaturação

**Padrões importantes:**
- **Oscilações regulares:** Dinâmica normal da proteína
- **Tendência decrescente:** Possível perda de estrutura
- **Mudanças abruptas:** Eventos de dobramento/desdobramento

## Scripts de Análise Avançada

### Análise de Convergência

```python
def analyze_convergence(data, window_size=100):
    """Analisa convergência de propriedades estruturais"""
    rolling_mean = np.convolve(data, np.ones(window_size)/window_size, mode='valid')
    slope, _, r_value, _, _ = stats.linregress(range(len(rolling_mean)), rolling_mean)
    
    convergence_score = abs(slope) / np.mean(data)
    
    if convergence_score < 0.001:
        return "Excelente convergência"
    elif convergence_score < 0.01:
        return "Boa convergência"
    else:
        return "Convergência questionável"
```

### Análise de Correlação

```python
def correlation_analysis(rmsd, gyration, hbonds):
    """Analisa correlações entre propriedades"""
    correlations = {
        'RMSD vs Gyration': np.corrcoef(rmsd, gyration)[0,1],
        'RMSD vs H-bonds': np.corrcoef(rmsd, hbonds)[0,1],
        'Gyration vs H-bonds': np.corrcoef(gyration, hbonds)[0,1]
    }
    
    for pair, corr in correlations.items():
        if abs(corr) > 0.7:
            print(f"{pair}: Forte correlação ({corr:.3f})")
        elif abs(corr) > 0.3:
            print(f"{pair}: Correlação moderada ({corr:.3f})")
        else:
            print(f"{pair}: Correlação fraca ({corr:.3f})")
```

## Checklist de Qualidade da Simulação

### ✅ Critérios de uma Simulação Bem-Sucedida

- [ ] **RMSD convergiu** para um valor estável (<0.3 nm para lisozima)
- [ ] **Raio de giração** apresenta baixa variação (<5%)
- [ ] **Ligações de hidrogênio** mantém-se estáveis (100-130 para lisozima)
- [ ] **RMSF** mostra padrão coerente (loops flexíveis, estruturas secundárias rígidas)
- [ ] **Energia** permanece estável durante a simulação produtiva
- [ ] **Temperatura e pressão** bem equilibradas

### ⚠️ Sinais de Alerta

- [ ] RMSD crescente continuamente
- [ ] Grandes flutuações no raio de giração (>10%)
- [ ] Perda significativa de ligações de hidrogênio (>20%)
- [ ] RMSF uniformemente alto (>0.4 nm para todos os resíduos)
- [ ] Energia instável ou crescente

## Exportação de Resultados

### Gerar Relatório Automatizado

```python
def generate_report(rmsd, rmsf, gyration, hbonds):
    """Gera relatório automatizado da simulação"""
    
    report = f"""
    RELATÓRIO DE SIMULAÇÃO DE DINÂMICA MOLECULAR
    ===========================================
    
    RESUMO ESTATÍSTICO:
    - RMSD médio: {np.mean(rmsd):.3f} ± {np.std(rmsd):.3f} nm
    - Raio de giração: {np.mean(gyration):.3f} ± {np.std(gyration):.3f} nm
    - Ligações H médias: {np.mean(hbonds):.1f} ± {np.std(hbonds):.1f}
    - Resíduos flexíveis: {len(rmsf[rmsf > 0.25])} de {len(rmsf)}
    
    AVALIAÇÃO DE QUALIDADE:
    {quality_assessment(rmsd, gyration, hbonds)}
    """
    
    return report
```

## Recursos Adicionais

### Visualização 3D com PyMOL

```python
# Script PyMOL para visualizar trajetória
pymol_script = '''
load md_0_1.pdb
load_traj md_0_1.xtc
color_by_rmsf rmsf.txt
show cartoon
set cartoon_transparency, 0.5
'''
```

### Análise com MDAnalysis

```python
import MDAnalysis as mda

# Carrega trajetória
u = mda.Universe('md_0_1.gro', 'md_0_1.xtc')

# Análises avançadas
protein = u.select_atoms('protein')
ca_atoms = u.select_atoms('protein and name CA')

# RMSF por resíduo
rmsf = rms.RMSF(ca_atoms).run()
```

Este guia fornece uma base sólida para interpretar e analisar os resultados de simulações de dinâmica molecular. Para dúvidas específicas ou análises mais avançadas, consulte a documentação do GROMACS e os artigos de referência listados no README principal.
