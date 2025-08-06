# Simulação de Dinâmica Molecular da Lisozima com GROMACS

Este repositório contém um tutorial completo e detalhado para a simulação de dinâmica molecular (MD) da lisozima em solução aquosa utilizando o GROMACS. A dinâmica molecular é uma técnica computacional indispensável para investigar a estrutura, a dinâmica e a energética de biomoléculas em resolução atômica. Este material abrange desde a preparação inicial do sistema até a simulação produtiva, com explicações e códigos detalhados para cada etapa.

---

Deixo um artigo de revisão sobre tutoriais de simulações em Dinâmica Molecular, como material de apoio.

```
Lemkul, J. A. Introductory Tutorials for Simulating Protein Dynamics with GROMACS. J. Phys. Chem. B 2024, 128 (39), 9418-9435. DOI: 10.1021/acs.jpcb.4c04901
```
Disponível para leitura em [PDF](https://pubs.acs.org/doi/pdf/10.1021/acs.jpcb.4c04901).

---

## Índice

- [Introdução](#introdução)
- [Requisitos](#requisitos)
- [Campos de Força (Force Fields)](#campos-de-força-force-fields)
- [Tipos de Arquivos GROMACS](#tipos-de-arquivos-gromacs)
- [Fluxo de Trabalho](#fluxo-de-trabalho)
- [Detalhamento das Etapas](#detalhamento-das-etapas)
  - [1. Preparação da Topologia da Proteína](#1-preparação-da-topologia-da-proteína)
  - [2. Definição da Caixa e Solvatação](#2-definição-da-caixa-e-solvatação)
  - [3. Adição de Íons](#3-adição-de-íons)
  - [4. Minimização de Energia](#4-minimização-de-energia)
  - [5. Equilibração NVT (Temperatura Constante)](#5-equilibração-nvt-temperatura-constante)
  - [6. Equilibração NPT (Pressão Constante)](#6-equilibração-npt-pressão-constante)
  - [7. Simulação Produtiva (MD)](#7-simulação-produtiva-md)
  - [8. Pós-processamento e Análise de Resultados](#8-pós-processamento-e-análise-de-resultados)
- [Como Executar o Tutorial](#como-executar-o-tutorial)
- [Possíveis Erros e Soluções](#possíveis-erros-e-soluções)
- [Referências](#referências)
- [Licença](#licença)

---

## Introdução

A dinâmica molecular permite estudar sistemas biomoleculares com uma resolução espacial e temporal que excede a maioria dos métodos experimentais. A lisozima é um sistema modelo clássico para esses estudos devido à sua estabilidade e estrutura bem caracterizada. Este tutorial guia o usuário através de um fluxo de trabalho fundamental para simular a lisozima em água, estabelecendo uma base sólida para projetos de simulação mais complexos.

---

## Requisitos

Antes de iniciar, certifique-se de ter os seguintes softwares instalados:

- **GROMACS** (versão recomendada: 2022 ou superior)
- **Python** (para scripts de análise e automação)
- **Visualizador Molecular:** `VMD` ou `PyMOL` para inspeção de trajetórias.
- **Arquivo PDB:** Estrutura da lisozima (ex: `1AKI.pdb` do Protein Data Bank).

Instale o GROMACS com:
```bash
sudo apt install gromacs
```
Ou, preferencialmente, via Conda para um melhor gerenciamento de ambientes:
```bash
conda install -c bioconda gromacs
```

---

## Campos de Força (Force Fields)

A precisão de uma simulação de MD depende fundamentalmente da qualidade do campo de força (force field) utilizado. Um campo de força é um conjunto de equações e parâmetros associados que descrevem a energia de uma dada configuração de átomos, aplicando os princípios da mecânica newtoniana para prever a evolução do sistema ao longo do tempo.

Existem muitos campos de força disponíveis, como **CHARMM**, **AMBER**, **GROMOS** e **OPLS**. Cada um possui parametrizações e pressupostos específicos. É crucial que as configurações da simulação (especificadas no arquivo `.mdp`) sejam consistentes com o campo de força escolhido. Isso inclui:

- **Modelo de Água:** Cada campo de força é calibrado com um modelo de água específico. Por exemplo, o campo de força CHARMM36 utiliza uma versão modificada do TIP3P. Usar o modelo de água incorreto pode levar a resultados imprecisos.
- **Parâmetros de Interação:** As configurações para o tratamento de interações não-ligadas (eletrostática e van der Waals), como os raios de corte (`cutoff`), devem ser consideradas como parte integrante do campo de força e não devem ser alteradas sem uma justificativa robusta.

Ao usar o `gmx pdb2gmx`, você selecionará um dos campos de força disponíveis para o GROMACS, e todas as etapas subsequentes devem respeitar as convenções desse campo de força.

---

## Tipos de Arquivos GROMACS

O GROMACS utiliza uma variedade de tipos de arquivos com extensões específicas. Abaixo estão os mais comuns encontrados neste tutorial:

- **`.gro`:** um arquivo de coordenadas de formato fixo com coordenadas dadas em unidades de nm.
- **`.pdb`:** um arquivo de coordenadas de formato fixo usado pelo Protein Databank com coordenadas em unidades de Å.
- **`.top`:** uma topologia do sistema, definindo o conteúdo completo de um sistema.
- **`.itp`:** uma topologia "incluída", definindo um tipo de molécula específico, parâmetros auxiliares ou outras diretivas topológicas.
- **`.mdp`:** arquivo de "parâmetros de dinâmica molecular" que especifica todas as configurações relevantes para realizar um cálculo ou simulação.
- **`.tpr`:** um arquivo de entrada de execução binário que combina coordenadas, topologia, todos os parâmetros do campo de força associados e todas as configurações de entrada definidas no arquivo .mdp.
- **`.edr`:** um arquivo binário contendo dados de energia do cálculo ou simulação.
- **`.xtc`:** um arquivo de trajetória binário em formato compactado contendo informações de tempo, vetor da caixa e coordenadas.
- **`.trr`:** um arquivo de trajetória de alta precisão contendo informações de tempo, vetor da caixa, coordenadas, velocidade e força.

---

## Fluxo de Trabalho

1.  **Preparação da Topologia:** Gerar uma topologia compatível com o campo de força a partir da estrutura PDB.
2.  **Definição da Caixa e Solvatação:** Criar uma caixa de simulação periódica e preenchê-la com moléculas de solvente (água).
3.  **Adição de Íons:** Adicionar íons para neutralizar a carga líquida do sistema e atingir uma concentração iônica desejada.
4.  **Minimização de Energia:** Remover contatos estéricos desfavoráveis e relaxar a geometria do sistema.
5.  **Equilibração NVT:** Estabilizar a temperatura do sistema, mantendo volume e número de partículas constantes.
6.  **Equilibração NPT:** Estabilizar a pressão e a densidade do sistema, mantendo temperatura e número de partículas constantes.
7.  **Simulação Produtiva:** Executar a simulação para coletar os dados, sem o uso de restrições.
8.  **Pós-processamento e Análise:** Extrair e analisar as propriedades estruturais e dinâmicas da proteína.

---

## Detalhamento das Etapas

### 1. Preparação da Topologia da Proteína

O primeiro passo é usar `gmx pdb2gmx` para ler o arquivo PDB (aqui chamado `model.pdb`), gerar uma topologia (`topol.top`) e um arquivo de coordenadas (`.gro`) compatível com o GROMACS.

```bash
# O GROMACS solicitará que você selecione um campo de força e um modelo de água da lista.
gmx pdb2gmx -f model.pdb -o model_processed.gro -water spce -ignh
```
- **Seleção Interativa:**
  1.  Escolha um campo de força da lista (ex: `CHARMM36-jul2022`).
  2.  Escolha o modelo de água correspondente (ex: `TIP3P`).

### 2. Definição da Caixa e Solvatação

Definimos uma caixa de simulação periódica e a preenchemos com moléculas de água.

```bash
# Define uma caixa cúbica a 1.0 nm da superfície da proteína.
gmx editconf -f model_processed.gro -o model_newbox.gro -c -d 1.0 -bt cubic
```bash
# Preenche a caixa com moléculas de água.
gmx solvate -cp model_newbox.gro -cs spc216.gro -o model_solv.gro -p topol.top
```

### 3. Adição de Íons

Neutralizamos a carga do sistema e adicionamos íons para simular uma concentração fisiológica.

```bash
# Prepara o sistema para a adição de íons.
gmx grompp -f ions.mdp -c model_solv.gro -p topol.top -o ions.tpr -maxwarn 5
```bash
# Adiciona íons Na+ e Cl- a uma concentração de 0.15 M e neutraliza a carga do sistema.
gmx genion -s ions.tpr -o model_solv_ions.gro -p topol.top -neutral -conc 0.15 -pname NA -nname CL
```
- **Seleção Interativa:** Quando solicitado, selecione o grupo de solvente a ser substituído pelos íons (geralmente `13 SOL`).

### 4. Minimização de Energia

Remove clivagens estéricas e relaxa a geometria do sistema.

```bash
# Prepara o sistema para a minimização.
gmx grompp -f minim.mdp -c model_solv_ions.gro -p topol.top -o em.tpr -maxwarn 5

# Executa a minimização de energia.
gmx mdrun -v -deffnm em -s em.tpr
```

### 5. Equilibração NVT (Temperatura Constante)

Estabiliza a temperatura do sistema com restrições de posição na proteína.

```bash
# Prepara o sistema para o equilíbrio NVT.
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro -maxwarn 5

# Executa o equilíbrio NVT.
gmx mdrun -deffnm nvt -v -s nvt.tpr
```

### 6. Equilibração NPT (Pressão Constante)

Estabiliza a pressão e a densidade do sistema, mantendo as restrições.

```bash
# Prepara o sistema para o equilíbrio NPT.
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -r em.gro -maxwarn 5

# Executa o equilíbrio NPT.
gmx mdrun -deffnm npt -v -s npt.tpr
```

### 7. Simulação Produtiva (MD)

Executa a simulação principal sem restrições para coletar dados.

```bash
# Prepara o sistema para a simulação de produção.
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr -maxwarn 5

# Executa a simulação produtiva.
gmx mdrun -deffnm md_0_1 -v
```

### 8. Pós-processamento e Análise de Resultados

Analisa a trajetória gerada (`md_0_1.xtc`).

#### 8.1 Cálculo do RMSF (Root-Mean-Square Fluctuation)

```bash
# O GROMACS solicitará a seleção do grupo para o cálculo.
gmx rmsf -s md_0_1.tpr -f md_0_1.xtc -o rmsf.xvg -res
```
- **Seleção Interativa:** Escolha `3 C-alpha`.

#### 8.2 Cálculo do RMSD (Root-Mean-Square Deviation)

```bash
# O GROMACS solicitará a seleção do grupo para o ajuste e para o cálculo.
gmx rms -s md_0_1.tpr -f md_0_1.xtc -o rmsd.xvg -tu ns
```
- **Seleção Interativa:** Escolha `4 Backbone` para o ajuste e `4 Backbone` novamente para o cálculo.

#### 8.3 Cálculo do Raio de Giração (Radius of Gyration)

```bash
# O GROMACS solicitará a seleção do grupo para o cálculo.
gmx gyrate -s md_0_1.tpr -f md_0_1.xtc -o giracao.xvg
```
- **Seleção Interativa:** Escolha `1 Protein`.

#### 8.4 Cálculo do Número de Ligações de Hidrogênio

```bash
# O GROMACS solicitará a seleção de dois grupos para o cálculo.
gmx hbond -s md_0_1.tpr -f md_0_1.xtc -num hbond_intra.xvg
```
- **Seleção Interativa:** Escolha `1 Protein` para o primeiro grupo e `1 Protein` novamente para o segundo.

#### 8.5 Visualização dos Resultados com Python

Para uma análise mais completa, você pode usar Python para gerar gráficos dos arquivos `.xvg` gerados pelo GROMACS:

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Função para ler arquivos .xvg do GROMACS
def read_xvg(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Remove linhas de comentário
    data_lines = [line for line in lines if not line.startswith('#') and not line.startswith('@')]
    
    # Converte para numpy array
    data = np.array([line.split() for line in data_lines], dtype=float)
    return data

# Configurações gerais dos gráficos
plt.style.use('seaborn-v0_8')
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
fig.suptitle('Análise de Dinâmica Molecular - Lisozima', fontsize=16, fontweight='bold')

# 1. RMSD (Root-Mean-Square Deviation)
rmsd_data = read_xvg('rmsd.xvg')
axes[0,0].plot(rmsd_data[:,0], rmsd_data[:,1], 'b-', linewidth=1.5)
axes[0,0].set_xlabel('Tempo (ns)')
axes[0,0].set_ylabel('RMSD (nm)')
axes[0,0].set_title('RMSD do Backbone da Proteína')
axes[0,0].grid(True, alpha=0.3)

# 2. RMSF (Root-Mean-Square Fluctuation)
rmsf_data = read_xvg('rmsf.xvg')
axes[0,1].plot(rmsf_data[:,0], rmsf_data[:,1], 'r-', linewidth=1.5)
axes[0,1].set_xlabel('Número do Resíduo')
axes[0,1].set_ylabel('RMSF (nm)')
axes[0,1].set_title('Flutuação por Resíduo (C-alpha)')
axes[0,1].grid(True, alpha=0.3)

# 3. Raio de Giração
gyration_data = read_xvg('giracao.xvg')
axes[1,0].plot(gyration_data[:,0], gyration_data[:,1], 'g-', linewidth=1.5)
axes[1,0].set_xlabel('Tempo (ns)')
axes[1,0].set_ylabel('Raio de Giração (nm)')
axes[1,0].set_title('Raio de Giração da Proteína')
axes[1,0].grid(True, alpha=0.3)

# 4. Ligações de Hidrogênio
hbond_data = read_xvg('hbond_intra.xvg')
axes[1,1].plot(hbond_data[:,0], hbond_data[:,1], 'm-', linewidth=1.5)
axes[1,1].set_xlabel('Tempo (ns)')
axes[1,1].set_ylabel('Número de Ligações H')
axes[1,1].set_title('Ligações de Hidrogênio Intramoleculares')
axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('analise_md_completa.png', dpi=300, bbox_inches='tight')
plt.show()

# Análise estatística básica
print("=== RESUMO ESTATÍSTICO ===")
print(f"RMSD médio: {np.mean(rmsd_data[:,1]):.3f} ± {np.std(rmsd_data[:,1]):.3f} nm")
print(f"Raio de giração médio: {np.mean(gyration_data[:,1]):.3f} ± {np.std(gyration_data[:,1]):.3f} nm")
print(f"Ligações H médias: {np.mean(hbond_data[:,1]):.1f} ± {np.std(hbond_data[:,1]):.1f}")
```

#### 8.6 Interpretação dos Resultados

**RMSD (Root-Mean-Square Deviation):**

- Valores típicos para proteínas estáveis: 0.1-0.3 nm
- Tendência crescente indica possível desnaturação
- Plateau indica estabilização da estrutura

![Exemplo RMSD](https://via.placeholder.com/600x400/4CAF50/FFFFFF?text=RMSD+vs+Tempo)

**RMSF (Root-Mean-Square Fluctuation):**

- Loops e terminais: alta flexibilidade (>0.3 nm)
- Folhas-β e α-hélices: baixa flexibilidade (<0.2 nm)
- Picos indicam regiões móveis importantes

![Exemplo RMSF](https://via.placeholder.com/600x400/2196F3/FFFFFF?text=RMSF+por+Resíduo)

**Raio de Giração:**

- Proteína compacta: ~1.4-1.6 nm para lisozima
- Variações pequenas (<5%) indicam estabilidade
- Aumentos significativos sugerem desdobramento

![Exemplo Raio de Giração](https://via.placeholder.com/600x400/FF9800/FFFFFF?text=Raio+de+Giração+vs+Tempo)

**Ligações de Hidrogênio:**

- Lisozima: ~100-130 ligações H intramoleculares
- Flutuações normais: ±10-15 ligações
- Perdas significativas indicam instabilidade estrutural

![Exemplo Ligações H](https://via.placeholder.com/600x400/9C27B0/FFFFFF?text=Ligações+H+vs+Tempo)

---

## Como Executar o Tutorial

1.  **Clone o repositório:**
    ```bash
    git clone [https://github.com/madsondeluna/aula_dinamica_molecular.git](https://github.com/madsondeluna/aula_dinamica_molecular.git)
    ```
2.  **Acesse o diretório:**
    ```bash
    cd aula_dinamica_molecular
    ```
3.  **Execute os comandos** na sequência apresentada na seção [Detalhamento das Etapas](#detalhamento-das-etapas).

---

## Possíveis Erros e Soluções

1.  **Erro `Fatal error: number of coordinates in coordinate file does not match topology`**
    -   **Causa**: O arquivo de coordenadas (`.gro`) e o de topologia (`.top`) estão dessincronizados. Isso geralmente acontece se um passo que adiciona/remove moléculas (como `gmx solvate` ou `gmx genion`) foi executado sem a flag `-p topol.top`.
    -   **Solução**: Refaça o passo problemático garantindo que a flag `-p` seja usada para atualizar a topologia.

2.  **Mensagem `Water molecule cannot be settled` ou instabilidade na simulação (LINCS warnings)**
    -   **Causa**: Uma geometria inicial muito ruim ou um passo de tempo (`dt`) muito grande no arquivo `.mdp`.
    -   **Solução**: Certifique-se de que a minimização de energia foi bem-sucedida. Se o erro persistir, reduza o `dt` no seu arquivo `.mdp` (ex: de 0.002 para 0.001).

3.  **Erro de `pdb2gmx` sobre nomes de átomos ou resíduos não reconhecidos**
    -   **Causa**: O arquivo PDB contém nomes de átomos ou resíduos não padrão.
    -   **Solução**: Use a flag `-ignh` para que o GROMACS reconstrua todos os hidrogênios. Se o erro for com átomos pesados, pode ser necessário editar o arquivo PDB manualmente.

4.  **Visualização "quebrada" da molécula no VMD/PyMOL**
    -   **Causa**: Artefato das condições de contorno periódicas.
    -   **Solução**: Use `gmx trjconv` para processar a trajetória antes de visualizar:
        ```bash
        gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_whole.xtc -pbc whole
        gmx trjconv -s md_0_1.tpr -f md_whole.xtc -o md_final.xtc -pbc nojump -center
        ```
    - **Seleção Interativa:** Para ambos os comandos, selecione `1 Protein` para centrar e `0 System` para a saída.

---

## Referências

- Lemkul, J. A. *J. Phys. Chem. B* **2024**, *128*, 9418-9435.
- [GROMACS Manual](http://www.gromacs.org/Documentation)
- [MD Tutorials – Lysozyme in Water](https://www.mdtutorials.com/gmx/lysozyme/)

---

## Licença

Este projeto está licenciado sob a [MIT License](LICENSE).
