# Simulação de Dinâmica Molecular da Lisozima com GROMACS

Este repositório contém um tutorial completo e detalhado para a simulação de dinâmica molecular (MD) da lisozima em solução aquosa utilizando o GROMACS. A dinâmica molecular é uma técnica computacional indispensável para investigar a estrutura, a dinâmica e a energética de biomoléculas em resolução atômica. Este material abrange desde a preparação inicial do sistema até a simulação produtiva, com explicações e códigos detalhados para cada etapa.

---

Deixo um artigo de revisão sobre tutoriais de simulações em Dinâmica Molecular, como material de apoio.

```
Lemkul, J. A. Introductory Tutorials for Simulating Protein Dynamics with GROMACS. J. Phys. Chem. B 2024, 128 (39), 9418-9435. DOI: 10.1021/acs.jpcb.4c04901
```
Disponível para leitura em: [GitHub deste tutorial](https://github.com/madsondeluna/aula_dinamica_molecular/blob/main/lemkul-2024-introductory-tutorials-for-simulating-protein-dynamics-with-gromacs.pdf).*

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
- [Estrutura do Projeto](#estrutura-do-projeto)
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

- **.gro:** um arquivo de coordenadas de formato fixo com coordenadas dadas em unidades de nm.
- **.pdb:** um arquivo de coordenadas de formato fixo usado pelo Protein Databank com coordenadas em unidades de Å.
- **.top:** uma topologia do sistema, definindo o conteúdo completo de um sistema.
- **.itp:** uma topologia "incluída", definindo um tipo de molécula específico, parâmetros auxiliares ou outras diretivas topológicas.
- **.mdp:** arquivo de "parâmetros de dinâmica molecular" que especifica todas as configurações relevantes para realizar um cálculo ou simulação.
- **.tpr:** um arquivo de entrada de execução binário que combina coordenadas, topologia, todos os parâmetros do campo de força associados e todas as configurações de entrada definidas no arquivo .mdp.
- **.edr:** um arquivo binário contendo dados de energia do cálculo ou simulação.
- **.xtc:** um arquivo de trajetória binário em formato compactado contendo informações de tempo, vetor da caixa e coordenadas.
- **.trr:** um arquivo de trajetória de alta precisão contendo informações de tempo, vetor da caixa, coordenadas, velocidade e força.

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

O primeiro passo é usar `gmx pdb2gmx` para ler o arquivo PDB, gerar uma topologia (`topol.top`) e um arquivo de coordenadas (`.gro`) compatível com o GROMACS. A topologia descreve as interações moleculares (ligações, ângulos, etc.) com base em um campo de força selecionado. O programa também adiciona átomos de hidrogênio que faltam na estrutura cristalina.

```bash
# Selecione interativamente um campo de força (ex: CHARMM36) e um modelo de água correspondente (ex: TIP3P).
gmx pdb2gmx -f 1AKI.pdb -o 1AKI_processed.gro -p topol.top -ignh -water tip3p
```
- `-f`: Arquivo de entrada PDB.
- `-o`: Arquivo de saída de coordenadas no formato GROMACS.
- `-p`: Arquivo de saída da topologia do sistema.
- `-water`: Especifica o modelo de água a ser usado, que deve ser compatível com o campo de força escolhido.
- `-ignh`: Ignora os hidrogênios do PDB e os reconstrói, útil para evitar erros de nomenclatura.

Este comando também cria o arquivo `posre.itp`, que contém definições para restrições de posição usadas durante a equilibração.

### 2. Definição da Caixa e Solvatação

Primeiro, definimos uma caixa de simulação periódica em torno da proteína com `gmx editconf`. Em seguida, preenchemos o volume vazio dessa caixa com moléculas de água usando `gmx solvate`.

```bash
# Define uma caixa (ex: dodecaedro rômbico) a 1.2 nm da superfície da proteína.
gmx editconf -f 1AKI_processed.gro -o 1AKI_box.gro -c -d 1.2 -bt dodecahedron
```
- `-c`: Centra a proteína na caixa.
- `-d`: Define a distância mínima entre a proteína e a borda da caixa. 1.0-1.2 nm é um valor comum.
- `-bt`: Define o tipo de caixa. O `dodecahedron` é mais eficiente em termos de volume que o `cubic` para solutos globulares, economizando custo computacional.

```bash
# Preenche a caixa com moléculas de água (spc216.gro é um arquivo genérico de solvente).
gmx solvate -cp 1AKI_box.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```
- `-cp`: Coordenadas do soluto (a caixa com a proteína).
- `-cs`: Coordenadas do solvente.
- `-p`: **Importante:** Atualiza automaticamente o arquivo `topol.top` com o número de moléculas de água adicionadas.

### 3. Adição de Íons

Para neutralizar a carga total do sistema (um requisito para o algoritmo PME de cálculo de eletrostática) e simular condições iônicas fisiológicas, adicionamos íons.

```bash
# Primeiro, crie um arquivo .tpr binário com o gmx grompp. Este arquivo contém toda a informação do sistema.
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr -maxwarn 1
```
- `ions.mdp`: Arquivo de parâmetros mínimo, pode estar vazio ou conter apenas `integrator = steep`.
- O `grompp` combina coordenadas, topologia e parâmetros em um único arquivo de entrada binário `.tpr`.
- `-maxwarn 1`: Permite ignorar avisos comuns nesta etapa.

```bash
# Use gmx genion para substituir moléculas de água por íons.
# O programa solicitará a seleção do grupo a ser substituído (geralmente "13 SOL" para a água).
echo "SOL" | gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```
- `-s`: Arquivo de entrada `.tpr`.
- `-pname` e `-nname`: Nomes do cátion e ânion, respectivamente.
- `-neutral`: Adiciona automaticamente o número de contra-íons necessários para zerar a carga líquida do sistema.
- `echo "SOL"`: Automatiza a seleção do grupo de solvente a ser substituído pelos íons.

### 4. Minimização de Energia

A minimização remove clivagens estéricas e geometrias desfavoráveis geradas durante a preparação do sistema, movendo os átomos para um mínimo de energia local. Isso garante um ponto de partida estável para a dinâmica.

```bash
# Prepara o sistema para a minimização.
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

# Executa a minimização de energia.
gmx mdrun -v -deffnm em
```
- `-deffnm em`: Uma abreviação conveniente que define `em` como o nome base para todos os arquivos de entrada e saída (`em.tpr`, `em.log`, `em.gro`, etc.).

### 5. Equilibração NVT (Temperatura Constante)

A primeira fase de equilibração estabiliza a temperatura do sistema. Restrições de posição são aplicadas aos átomos pesados da proteína para permitir que as moléculas de solvente se acomodem ao redor dela sem perturbar sua estrutura.

```bash
# Prepara o sistema para o equilíbrio NVT.
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

# Executa o equilíbrio NVT.
gmx mdrun -deffnm nvt
```
- `-r em.gro`: Especifica as coordenadas de referência para as restrições de posição. As restrições são ativadas no arquivo `nvt.mdp` pela linha `define = -DPOSRES`, que por sua vez inclui o arquivo `posre.itp` na topologia.

### 6. Equilibração NPT (Pressão Constante)

A segunda fase de equilibração estabiliza a pressão, ajustando o volume da caixa para atingir a densidade correta do sistema. As restrições de posição na proteína são mantidas.

```bash
# Prepara o sistema para o equilíbrio NPT, usando o estado final do NVT como ponto de partida.
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

# Executa o equilíbrio NPT.
gmx mdrun -deffnm npt
```
- `-c nvt.gro`: Usa as coordenadas finais da etapa NVT.
- `-t nvt.cpt`: Usa o arquivo de checkpoint do NVT para uma continuação exata, preservando velocidades e estado do termostato/barostato.

### 7. Simulação Produtiva (MD)

Após a equilibração, as restrições são removidas e a simulação é executada para coletar os dados que serão analisados. A duração depende do fenômeno de interesse, variando de nanossegundos (ns) a microssegundos (µs).

```bash
# Prepara o sistema para a simulação de produção.
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

# Executa a simulação produtiva.
gmx mdrun -deffnm md
```
- O arquivo `md.mdp` é semelhante ao `npt.mdp`, mas sem a linha `define = -DPOSRES`.

### 8. Pós-processamento e Análise de Resultados

Antes da análise, é crucial corrigir os efeitos de borda periódica na trajetória para visualização e cálculo corretos. Use `gmx trjconv` para garantir que a proteína não apareça "quebrada" ou "pulando" na caixa.

#### 8.1 Cálculo do RMSF (Root-Mean-Square Fluctuation)

Mede a flutuação de cada resíduo em torno de sua posição média, indicando as regiões mais flexíveis da proteína.

```bash
# Selecione o grupo "C-alpha" (geralmente 3) para o cálculo.
echo "C-alpha" | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res
```
- `-res`: Calcula a flutuação por resíduo, não por átomo.

#### 8.2 Cálculo do RMSD (Root-Mean-Square Deviation)

Mede o desvio estrutural global da proteína ao longo do tempo em relação a uma estrutura de referência (ex: a estrutura inicial). É um indicador da estabilidade conformacional.

```bash
# Selecione "Backbone" (esqueleto) para o ajuste e "Backbone" para o cálculo.
echo "Backbone Backbone" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns
```
- `-tu ns`: Exibe o tempo no eixo X em nanossegundos.

#### 8.3 Cálculo do Raio de Giração (Radius of Gyration)

Mede a compactação da proteína. Aumentos no raio de giração podem indicar um processo de desdobramento (unfolding).

```bash
# Selecione "Protein" (geralmente 1) para o cálculo.
echo "Protein" | gmx gyrate -s md.tpr -f md.xtc -o giracao.xvg
```

#### 8.4 Cálculo do Número de Ligações de Hidrogênio

Analisa as ligações de hidrogênio, que são fundamentais para a manutenção da estrutura secundária e terciária da proteína.

```bash
# Para ligações de H intramoleculares, selecione "Protein" duas vezes.
echo "Protein Protein" | gmx hbond -s md.tpr -f md.xtc -num hbond_intra.xvg
```

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
    -   **Causa**: Uma geometria inicial muito ruim (após adicionar íons, por exemplo) ou um passo de tempo (`dt`) muito grande no arquivo `.mdp`.
    -   **Solução**: Certifique-se de que a minimização de energia foi bem-sucedida. Se o erro persistir, reduza o `dt` no seu arquivo `.mdp` (ex: de 0.002 para 0.001) e tente novamente.

3.  **Erro de `pdb2gmx` sobre nomes de átomos ou resíduos não reconhecidos**
    -   **Causa**: O arquivo PDB contém nomes de átomos ou resíduos não padrão, ou hidrogênios com nomenclatura incorreta.
    -   **Solução**: Use a flag `-ignh` para que o GROMACS reconstrua todos os hidrogênios. Se o erro for com átomos pesados, pode ser necessário editar o arquivo PDB manualmente para corrigir os nomes.

4.  **Visualização "quebrada" da molécula no VMD/PyMOL**
    -   **Causa**: Isso não é um erro, mas um artefato das condições de contorno periódicas. A molécula pode cruzar a fronteira da caixa.
    -   **Solução**: Use `gmx trjconv` para processar a trajetória antes de visualizar. Centralize a proteína e torne as moléculas inteiras:
        ```bash
        # Selecione "Protein" para centrar e "System" para a saída.
        echo "Protein System" | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole
        echo "Protein System" | gmx trjconv -s md.tpr -f md_whole.xtc -o md_final.xtc -pbc nojump -center
        ```

---

## Referências

- Lemkul, J. A. *J. Phys. Chem. B* **2024**, *128*, 9418-9435.
- [GROMACS Manual](http://www.gromacs.org/Documentation)
- [MD Tutorials – Lysozyme in Water](https://www.mdtutorials.com/gmx/lysozyme/)

---

## Agradecimentos

- Um agradecimento especial ao projeto *Making-It-Rain* de Pablo R. Arantes (@pablitoarantes), Marcelo D. Polêto (@mdpoleto), Conrado Pedebos (@ConradoPedebos) e Rodrigo Ligabue-Braun (@ligabue_braun), que permitiu a execução dos protocolos na nuvem com o Colab.

---

## Licença

Este projeto está licenciado sob a [MIT License](LICENSE).
