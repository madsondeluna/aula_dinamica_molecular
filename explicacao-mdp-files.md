# Configuração de Arquivos `.mdp` no GROMACS

Os arquivos `.mdp` (Molecular Dynamics Parameter) contêm todos os parâmetros necessários para execução de uma simulação no GROMACS. Cada etapa do protocolo (minimização, NVT, NPT e MD produtiva) tem seu próprio arquivo `.mdp` com configurações específicas.

> **Importante:** Os parâmetros de não-ligados (cutoffs, vdW modifier) fazem parte do **campo de força** e não devem ser alterados sem justificativa robusta. As configurações abaixo são específicas para o **CHARMM36m**.

---

## Parâmetros por Categoria

### 1. Integrador e Controle de Tempo

| Parâmetro | Descrição | Valores comuns |
|-----------|-----------|----------------|
| `integrator` | Algoritmo de integração | `md` (leap-frog), `steep` (minimização) |
| `dt` | Passo de tempo | `0.002` ps (2 fs) para MD com restrições em H |
| `nsteps` | Número de passos | `50000` (100 ps), `250000` (500 ps), `25000000` (50 ns) |

### 2. Controle de Saída

| Parâmetro | Descrição |
|-----------|-----------|
| `nstxout-compressed` | Frequência de escrita do `.xtc` (trajetória comprimida) |
| `nstenergy` | Frequência de escrita do `.edr` (energia) |
| `nstlog` | Frequência de atualização do `.log` |

### 3. Interações Não-Ligadas (CHARMM36m)

```
cutoff-scheme   = Verlet        ; lista de vizinhos com buffer automatico
vdwtype         = cutoff
vdw-modifier    = force-switch  ; chaveamento suave: essencial para CHARMM36m
rvdw-switch     = 1.0           ; inicio do chaveamento (nm)
rvdw            = 1.2           ; cutoff vdW (nm)
rlist           = 1.4           ; cutoff lista de vizinhos (nm)
dispcorr        = no            ; CHARMM36m: sem correcao de dispersao
coulombtype     = PME           ; Particle Mesh Ewald para eletrostatica
rcoulomb        = 1.2           ; cutoff eletrostatico (nm)
fourierspacing  = 0.16          ; espacamento da grade FFT
```

> `vdw-modifier = force-switch` é a principal diferença entre CHARMM36m e OPLS/AMBER.
> Com OPLS/AMBER, usa-se `dispcorr = EnerPres` e `rvdw = 1.0`. **Não misture configurações de campos de força diferentes.**

### 4. Restrições de Ligação

```
constraint_algorithm = lincs    ; LINCS: mais estavel que SHAKE
constraints          = h-bonds  ; restringir apenas ligacoes com H (CHARMM36m)
lincs_iter           = 1
lincs_order          = 4
```

> CHARMM36m usa `constraints = h-bonds`. OPLS usava `all-bonds`. Usar `h-bonds` com CHARMM36m
> permite `dt = 0.002 ps` sem instabilidades.

### 5. Controle de Temperatura

```
tcoupl    = V-rescale    ; termostato de Bussi et al. (recomendado para MD)
tc-grps   = System       ; acoplar o sistema inteiro como um grupo
tau_t     = 1.0          ; constante de tempo (ps) — CHARMM36m usa 1.0, nao 0.1
ref_t     = 298          ; temperatura de referencia (K)
```

> CHARMM36m acopla o **System** como um único grupo (`tc-grps = System`), com `tau_t = 1.0 ps`.
> Configurações antigas do OPLS usavam `Protein Non-Protein` e `tau_t = 0.1 ps`.

### 6. Controle de Pressão

```
pcoupl          = C-rescale          ; barostato C-rescale (CHARMM36m, GROMACS 2021+)
pcoupltype      = isotropic
tau_p           = 5.0                ; constante de tempo (ps)
ref_p           = 1.0                ; pressao de referencia (bar)
compressibility = 4.5e-5             ; compressibilidade isot. da agua (bar^-1)
refcoord_scaling = com               ; escalar coordenadas de referencia (NVT/NPT com POSRES)
```

> O barostato **C-rescale** substituiu o Parrinello-Rahman como recomendação para CHARMM36m.
> O C-rescale é estocástico e mais estável durante a equilibração. Use `tau_p = 5.0 ps`.

### 7. Restrições de Posição

```
define = -DPOSRES    ; ativa as restricoes definidas em posre.itp
                     ; usar em NVT e NPT (remover na MD produtiva)
```

### 8. Geração de Velocidades

```
; Apenas no primeiro passo (NVT, continuation = no):
gen_vel  = yes
gen_temp = 298       ; deve ser igual a ref_t
gen_seed = -1        ; semente aleatoria

; Nas etapas seguintes (NPT, MD produtiva):
gen_vel  = no
continuation = yes
```

---

## Resumo dos Arquivos `.mdp` deste Repositório

| Arquivo | Etapa | nsteps | Destaque |
|---------|-------|--------|----------|
| `ions.mdp` | Geração de íons | 50000 | `coulombtype = cutoff` (apenas placeholder) |
| `minim.mdp` | Minimização | 50000 | `emtol = 500`, sem termostato/barostato |
| `nvt.mdp` | Equilibração NVT | 50000 (100 ps) | `define = -DPOSRES`, `gen_vel = yes`, sem barostato |
| `npt.mdp` | Equilibração NPT | 250000 (500 ps) | `define = -DPOSRES`, `pcoupl = C-rescale` |
| `md.mdp` | Produção | 25000000 (50 ns) | Sem restrições, `continuation = yes` |

---

## Comparação CHARMM36m vs. OPLS-AA (referência)

| Parâmetro | CHARMM36m | OPLS-AA |
|-----------|-----------|---------|
| `vdw-modifier` | `force-switch` | `none` (cut-off) |
| `rvdw` | 1.2 nm | 1.0 nm |
| `rcoulomb` | 1.2 nm | 1.0 nm |
| `dispcorr` | `no` | `EnerPres` |
| `constraints` | `h-bonds` | `all-bonds` |
| `tc-grps` | `System` | `Protein Non-Protein` |
| `tau_t` | 1.0 ps | 0.1 ps |
| `pcoupl` | `C-rescale` | `Parrinello-Rahman` |
| `tau_p` | 5.0 ps | 2.0 ps |
| `ref_t` | 298 K | 300 K |
