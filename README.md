# Simulação de Dinâmica Molecular da Lisozima com GROMACS

Este repositório contém um tutorial completo e detalhado para a simulação de dinâmica molecular (MD) da lisozima em solução aquosa utilizando o GROMACS. O material abrange desde a preparação inicial do sistema até a simulação produtiva, com explicações e códigos de cada etapa.

---

## Índice

- [Introdução](#introdução)
- [Fluxo de Trabalho](#fluxo-de-trabalho)
- [Detalhamento das Etapas](#detalhamento-das-etapas)
  - [1. Preparação do Sistema](#1-preparação-do-sistema)
  - [2. Criação da Caixa e Solvatação](#2-criação-da-caixa-e-solvatação)
  - [3. Adição de Íons](#3-adição-de-íons)
  - [4. Minimização de Energia](#4-minimização-de-energia)
  - [5. Equilíbrio NVT](#5-equilíbrio-nvt)
  - [6. Equilíbrio NPT](#6-equilíbrio-npt)
  - [7. Simulação](#7-simulação-produtiva)
  - [8. Pós-processamento e Análise de Resultados](#8-pós-processamento-e-análise-de-resultados)
- [Estrutura do Projeto](#estrutura-do-projeto)
- [Como Executar o Tutorial](#como-executar-o-tutorial)
- [Referências](#referências)
- [Licença](#licença)

---

## 1. Preparação do Sistema

```bash
gmx pdb2gmx -f 1AKI.pdb -o 1AKI_processed.gro -p topol.top -i posre.itp -water spce -ff oplsaa
```
- Converte o arquivo PDB para o formato GROMACS e gera os arquivos de topologia e restrições.

## 2. Criação da Caixa e Solvatação

```bash
gmx editconf -f 1AKI_processed.gro -o 1AKI_box.gro -c -d 1.0 -bt cubic
```
- Define uma caixa cúbica ao redor da proteína.

```bash
gmx solvate -cp 1AKI_box.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```
- Preenche a caixa com moléculas de água.

## 3. Adição de Íons

```bash
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
```
- Prepara o sistema para a adição de íons.

```bash
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```
- Substitui moléculas de água por íons para neutralizar o sistema.

## 4. Minimização de Energia

```bash
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```
- Prepara o sistema para a minimização de energia.

```bash
gmx mdrun -v -deffnm em
```
- Executa a minimização de energia.

## 5. Equilíbrio NVT

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```
- Prepara o sistema para o equilíbrio a temperatura constante.

```bash
gmx mdrun -deffnm nvt
```
- Executa o equilíbrio NVT.

## 6. Equilíbrio NPT

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
```
- Prepara o sistema para o equilíbrio a pressão constante.

```bash
gmx mdrun -deffnm npt
```
- Executa o equilíbrio NPT.

## 7. Simulação 

```bash
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
```
- Prepara o sistema para a simulação de produção.

```bash
gmx mdrun -deffnm md
```
- Executa a simulação produtiva.

## 8. Pós-processamento e Análise de Resultados

### 8.1 Cálculo do RMSD

```bash
gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res
```
- Mede a flutuação média dos resíduos da proteína ao longo da simulação.
### 8.2 Cálculo do RMSF

```bash
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns
```
- Mede o desvio médio quadrático das coordenadas da proteína.

### 8.3 Cálculo do Raio de Giração

```bash
gmx gyrate -s md.tpr -f md.xtc -o giracao.xvg
```
- Mede a compactação estrutural da proteína.

### 8.4 Cálculo do Número de Ligações de Hidrogênio

```bash
gmx hbond -s md.tpr -f md.xtc -num hbond.xvg
```
- Analisa a quantidade de ligações de hidrogênio formadas.

---

## Estrutura do Projeto

```
.
├── 1AKI.pdb
├── mdp/
│   ├── ions.mdp
│   ├── minim.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   ├── md.mdp
├── scripts/
│   ├── preparar.sh
│   ├── rodar_md.sh
│   ├── analisar.sh
├── README.md
```

---

## Como Executar o Tutorial

1. **Clone o repositório:**
   ```bash
   git clone https://github.com/seu_usuario/seu_repositorio.git
   ```
2. **Acesse o diretório:**
   ```bash
   cd seu_repositorio
   ```
3. **Execute os scripts disponíveis para automação do processo.**

---

## Referências

- [GROMACS Manual](http://www.gromacs.org/Documentation)
- [MD Tutorials – Lysozyme in Water](https://www.mdtutorials.com/gmx/lysozyme/)

---

## Licença

Este projeto está licenciado sob a [MIT License](LICENSE).

