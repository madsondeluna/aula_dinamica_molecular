# Tutoriais Introdutórios de Dinâmica Molecular com GROMACS

Este repositório contém três tutoriais práticos para simulação de dinâmica molecular (MD) de proteínas utilizando o GROMACS, baseados no artigo:

```
Lemkul, J. A. Introductory Tutorials for Simulating Protein Dynamics with GROMACS.
J. Phys. Chem. B 2024, 128 (39), 9418-9435. DOI: 10.1021/acs.jpcb.4c04901
```

Os arquivos de input originais estão disponíveis em: https://github.com/Lemkul-Lab/gmx_tutorials_jpcb

---

## Índice

- [Visão Geral](#visão-geral)
- [Pré-requisitos](#pré-requisitos)
- [Campo de Força CHARMM36m](#campo-de-força-charmm36m)
- [Tipos de Arquivos GROMACS](#tipos-de-arquivos-gromacs)
- [Exercício 1 — Proteína em Água (Ubiquitina)](#exercício-1--proteína-em-água-ubiquitina)
- [Exercício 2 — Complexo Proteico em Água (InaD:NorpA)](#exercício-2--complexo-proteico-em-água-inadnorpa)
- [Exercício 3 — Umbrella Sampling (Chignolin)](#exercício-3--umbrella-sampling-chignolin)
- [Possíveis Erros e Soluções](#possíveis-erros-e-soluções)
- [Referências](#referências)

---

## Visão Geral

Os três exercícios cobrem fluxos de trabalho fundamentais em MD:

| Exercício | Sistema | PDB | Conteúdo |
|-----------|---------|-----|----------|
| 1 | Ubiquitina em água | 1UBQ | Preparação, simulação e análise de uma proteína globular |
| 2 | Complexo InaD:NorpA | 1IHJ | Complexo proteico com ligação dissulfeto intermolecular e grupos de cap |
| 3 | Chignolin (β-hairpin) | 1UAO | Umbrella sampling e cálculo de energia livre (PMF/WHAM) |

---

## Pré-requisitos

- **GROMACS** 2024 ou superior
- **Campo de Força CHARMM36m** (julho 2022) — download externo (ver abaixo)
- **Python 3** com `matplotlib`, `numpy` e `pandas`
- **VMD** ou **PyMOL** para visualização de trajetórias
- **PyMOL** também é necessário no Exercício 2 (para adicionar grupos de cap)

**Instalar GROMACS via conda (recomendado):**
```bash
conda install -c conda-forge gromacs
```
**Ou via apt:**
```bash
sudo apt-get install gromacs
```

---

## Campo de Força CHARMM36m

Todos os exercícios utilizam o **CHARMM36m** (versão julho de 2022), que **não está incluído** na instalação padrão do GROMACS e deve ser baixado separadamente:

```bash
wget https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz \
     -O charmm36-jul2022.ff.tgz
tar xzf charmm36-jul2022.ff.tgz
```

Coloque o diretório `charmm36-jul2022.ff` no seu diretório de trabalho (onde os comandos GROMACS serão executados) antes de iniciar qualquer exercício.

### Configurações obrigatórias do CHARMM36m

As seguintes configurações devem ser usadas **em todas as etapas** (minimização, NVT, NPT, MD produtiva):

```
vdw-modifier    = force-switch   ; chaveamento suave das forças de van der Waals
rvdw-switch     = 1.0            ; distância de início do chaveamento (nm)
rvdw            = 1.2            ; cutoff van der Waals (nm)
rcoulomb        = 1.2            ; cutoff eletrostático (nm)
rlist           = 1.4            ; cutoff da lista de vizinhos (nm)
dispcorr        = no             ; sem correção de dispersão de longo alcance
constraints     = h-bonds        ; apenas ligações com H são restritas
```

---

## Tipos de Arquivos GROMACS

| Extensão | Descrição |
|----------|-----------|
| `.gro`   | Arquivo de coordenadas (unidades em nm) |
| `.pdb`   | Arquivo de coordenadas do Protein Data Bank (unidades em Å) |
| `.top`   | Topologia completa do sistema |
| `.itp`   | Topologia "incluída" — define um tipo de molécula ou parâmetros auxiliares |
| `.mdp`   | Parâmetros de dinâmica molecular — configurações de simulação |
| `.tpr`   | Arquivo de input binário (coordenadas + topologia + parâmetros + .mdp) |
| `.edr`   | Arquivo binário de energia |
| `.xtc`   | Arquivo de trajetória comprimido (tempo, caixa, coordenadas) |
| `.trr`   | Arquivo de trajetória de alta precisão (tempo, caixa, coordenadas, velocidades, forças) |
| `.cpt`   | Arquivo de checkpoint para reiniciar simulações |

---

## Exercício 1 — Proteína em Água (Ubiquitina)

**Proteína:** Ubiquitina — PDB: [1UBQ](https://www.rcsb.org/structure/1UBQ)

Este é o exercício mais detalhado, cobrindo todo o fluxo de trabalho: preparação, simulação e análise.

### Fluxo de Trabalho

```
1UBQ.pdb → topologia → caixa dodecaédrica → solvatação → adição de íons
         → minimização → NVT (100 ps) → NPT (500 ps) → MD (50+50 ns)
         → análise (RMSD, RMSF, ligações H, DSSP)
```

---

### Passo 1 — Preparação da Topologia

```bash
gmx pdb2gmx -f 1UBQ.pdb -o ubiquitin.gro -water tip3p
```

Quando solicitado:
1. Escolha **CHARMM36** (aparecerá como opção 1 se `charmm36-jul2022.ff` estiver no diretório de trabalho)
2. Escolha **TIP3P** (opção 1 — modelo padrão para CHARMM36)

Arquivos gerados:
- `ubiquitin.gro` — estrutura processada e protonada
- `topol.top` — topologia completa
- `posre.itp` — restrições de posição para a equilibração

---

### Passo 2 — Definição da Caixa e Solvatação

Definimos uma **caixa dodecaédrica rômbica** (mais eficiente que a cúbica — menor volume com a mesma distância periódica):

```bash
gmx editconf -f ubiquitin.gro -o ubiquitin_box.gro -c -d 1.2 -bt dodecahedron
```

Solvatação com água TIP3P:

```bash
gmx solvate -cp ubiquitin_box.gro -cs spc216.gro -o ubiquitin_solv.gro -p topol.top
```

> **Nota:** O arquivo `spc216.gro` (caixa de água pré-equilibrada) é utilizado para qualquer modelo de água de três pontos — o sistema converge para as propriedades do TIP3P durante a minimização e equilibração.

---

### Passo 3 — Adição de Íons

A ubiquitina possui carga líquida nula. Adicionamos 0.1 M de NaCl para simular condições fisiológicas:

```bash
gmx grompp -f ions.mdp -c ubiquitin_solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o ubiquitin_solv_ions.gro -p topol.top \
             -pname NA -nname CL -conc 0.1
```

---

### Passo 4 — Minimização de Energia

```bash
gmx grompp -f minim.mdp -c ubiquitin_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

Para extrair e verificar a energia potencial:

```bash
echo -e "11\n0" | gmx energy -f em.edr -o potential.xvg
```

> Critério de convergência: força máxima < 500 kJ/mol/nm (`emtol = 500.0`).

---

### Passo 5 — Equilibração NVT (100 ps)

Estabiliza a temperatura a 298 K com restrições de posição na proteína:

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -nb gpu
```

Para verificar a temperatura:

```bash
echo -e "16\n0" | gmx energy -f nvt.edr -o temperature.xvg
```

---

### Passo 6 — Equilibração NPT (500 ps)

Estabiliza a pressão a 1 bar (barostato C-rescale):

```bash
gmx grompp -f npt.mdp -c nvt.gro -r em.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -nb gpu
```

Para verificar a pressão:

```bash
echo -e "17\n0" | gmx energy -f npt.edr -o pressure.xvg
```

---

### Passo 7 — Simulação Produtiva (50 ns + extensão para 100 ns)

**Primeira metade (0–50 ns):**

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_50.tpr
gmx mdrun -deffnm md_0_50 -nb gpu
```

**Extensão para 100 ns (50–100 ns):**

```bash
gmx convert-tpr -s md_0_50.tpr -until 100000 -o md_50_100.tpr
gmx mdrun -s md_50_100.tpr -deffnm md_50_100 -cpi md_0_50.cpt -noappend -nb gpu
```

---

### Passo 8 — Pós-processamento e Análise

#### 8.1 Concatenação e reimagem da trajetória

```bash
# Concatenar as duas metades
gmx trjcat -f md_0_50.xtc md_50_100.xtc -o md_all.xtc

# Reimagem: colocar moléculas na imagem central, centrar a proteína
echo -e "1\n0" | gmx trjconv -s md_0_50.tpr -f md_all.xtc -o md_reimaged.xtc \
                 -pbc mol -ur compact -center

# Ajuste de rotação e translação global
echo -e "4\n1" | gmx trjconv -s md_0_50.tpr -f md_reimaged.xtc -o md_fit.xtc \
                 -fit rot+trans
```

#### 8.2 RMSD (Root-Mean-Square Deviation)

```bash
echo -e "4\n4" | gmx rms -s md_0_50.tpr -f md_fit.xtc -o rmsd.xvg -tu ns
```
- Selecione **Backbone** para ajuste e cálculo.
- Valores típicos para proteínas estáveis: 0.1–0.2 nm.

#### 8.3 RMSF (Root-Mean-Square Fluctuation)

```bash
echo "3" | gmx rmsf -s md_0_50.tpr -f md_fit.xtc -o rmsf.xvg -res
```
- Selecione **C-alpha**.
- Picos indicam regiões flexíveis (loops e terminais).

#### 8.4 Ligações de Hidrogênio

```bash
# Total de ligações H intramoleculares
echo -e "1\n1" | gmx hbond -s md_0_50.tpr -f md_fit.xtc -num hbond_total.xvg

# Ligações H da cadeia principal (requer MainChain+H)
echo -e "6\n6" | gmx hbond -s md_0_50.tpr -f md_fit.xtc -num hbond_backbone.xvg
```

> **Importante:** Para ligações H envolvendo a cadeia principal, use o grupo **MainChain+H** (não "Backbone"), pois os hidrogênios amídicos são necessários para o critério geométrico do cálculo.

#### 8.5 Estrutura Secundária (DSSP)

```bash
echo "1" | gmx dssp -s md_0_50.tpr -f md_fit.xtc -o dssp.dat
```

---

## Exercício 2 — Complexo Proteico em Água (InaD:NorpA)

**Sistema:** Complexo InaD:NorpA (domínio PDZ de InaD com peptídeo NorpA) — PDB: [1IHJ](https://www.rcsb.org/structure/1IHJ)

Este exercício expande o fluxo do Exercício 1 com os seguintes desafios adicionais:
- Extração de cadeias específicas de um arquivo PDB com múltiplas cadeias
- Modelagem de grupos de cap (ACE/NME) para termini incompletos com PyMOL
- Reconstrução de átomos de cadeias laterais ausentes
- Ligação dissulfeto **intermolecular** entre Cys31 (InaD) e Cys6 (NorpA)
- Uso de `-merge` e `-ter` no `pdb2gmx`
- Uso de `-neutral` no `genion` para neutralizar a carga líquida do complexo (−1)

**Comando de exemplo para pdb2gmx:**

```bash
gmx pdb2gmx -f 1ihj_chainAD_capped.pdb -o complex.gro -water tip3p -ter -merge all
```

> Consulte o artigo original para o procedimento detalhado de preparação da estrutura com PyMOL.

---

## Exercício 3 — Umbrella Sampling (Chignolin)

**Sistema:** Chignolin (β-hairpin de 10 resíduos) — PDB: [1UAO](https://www.rcsb.org/structure/1UAO)

Este exercício demonstra o uso de potencial de bias para calcular superfícies de energia livre (PMF) usando WHAM:

1. **Preparação** — Topologia e caixa dodecaédrica com `-box 7.0`
2. **Pulling** — Separação do Cα de Gly1 ao Cα de Gly10 (taxa: 0.005 nm/ps, k: 2000 kJ/mol/nm²)
3. **Janelas de umbrella** — 26 janelas de 0.5 a 3.0 nm (espaçamento: 0.1 nm)
4. **Análise WHAM**:

```bash
gmx wham -it tpr_files.dat -if pullf_files.dat -o pmf.xvg -hist hist.xvg -unit kJ
```

**Grupos de índice para a coordenada de reação:**

```
make_ndx > r 1 & a CA
make_ndx > name 2 Gly1_CA
make_ndx > r 10 & a CA
make_ndx > name 3 Gly10_CA
```

> Consulte o artigo original para o script Python de seleção de frames e as configurações completas do arquivo `.mdp` de pulling.

---

## Possíveis Erros e Soluções

1. **`Fatal error: number of coordinates in coordinate file does not match topology`**
   - **Causa:** Arquivo `.gro` e `.top` dessincronizados.
   - **Solução:** Refaça o passo que adicionou/removeu moléculas garantindo que a flag `-p topol.top` seja usada.

2. **`Water molecule cannot be settled` ou instabilidade (LINCS warnings)**
   - **Causa:** Geometria inicial ruim ou `dt` muito grande.
   - **Solução:** Verifique se a minimização convergiu. Reduza `dt` de 0.002 para 0.001 ps se necessário.

3. **`pdb2gmx` não reconhece nomes de átomos ou resíduos**
   - **Causa:** Nomenclatura não padrão no arquivo PDB.
   - **Solução:** Use a flag `-ignh` para reconstruir hidrogênios.

4. **Molécula "quebrada" na visualização (VMD/PyMOL)**
   - **Causa:** Artefato das condições de contorno periódicas.
   - **Solução:** Faça a reimagem da trajetória:
     ```bash
     echo -e "1\n0" | gmx trjconv -s md_0_50.tpr -f md_all.xtc -o md_reimaged.xtc \
                      -pbc mol -ur compact -center
     ```

5. **`-nb gpu` falha (sem GPU disponível)**
   - **Solução:** Remova a flag `-nb gpu` e use somente CPU:
     ```bash
     gmx mdrun -deffnm nvt -ntmpi 1
     ```

---

## Referências

- Lemkul, J. A. *J. Phys. Chem. B* **2024**, *128*, 9418–9435. https://doi.org/10.1021/acs.jpcb.4c04901
- Repositório de inputs do tutorial: https://github.com/Lemkul-Lab/gmx_tutorials_jpcb
- GROMACS Manual: https://manual.gromacs.org/
- CHARMM36m force field: http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs

---

## Licença

Este projeto está licenciado sob a [MIT License](LICENSE).
