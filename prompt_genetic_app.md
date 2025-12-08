# PROMPT: Genetic Analysis Desktop Application (Níveis 1-2)

## OBJETIVO GLOBAL

Gerar código Python completo e funcional para uma aplicação desktop (PyQt6) que permite:
1. Upload de dados genéticos brutos 23andMe
2. Matching automático contra base de dados GWAS Catalog
3. Cálculo de impact scores
4. Filtragem, busca e visualização de resultados

A aplicação deve ser:
- **COMPLETAMENTE FUNCIONAL**: Sem TODOs, placeholders ou pseudocode
- **OFFLINE**: Base de dados carregada localmente, sem dependências de internet
- **ROBUSTA**: Error handling completo, validação de inputs, logging
- **ESCALÁVEL**: Arquitetura permite extensão para Níveis 3-4 posteriormente
- **TESTÁVEL**: Código modular com testes unitários

---

## ESCOPO FUNCIONAL DETALHADO

### NÍVEL 1: MVP Básico

**Fluxo Utilizador:**
1. Utilizador seleciona ficheiro 23andMe raw data (formato .txt)
2. Aplicação faz parsing do ficheiro, extrai SNPs
3. Cada SNP é procurado na base de dados GWAS Catalog
4. Para cada SNP encontrado em GWAS, retorna associações trait
5. Exibe tabela com colunas: SNP ID | Gene | Trait | Genótipo Utilizador | Alelo Risco | P-value
6. Tabela é ordenável por qualquer coluna
7. Paginação para datasets grandes (50 resultados/página)
8. Status bar mostra número total de matches

**Validações Nível 1:**
- Ficheiro 23andMe: validar formato (RSID, chromosome, position, genotype)
- SNPs user: rejeitar linhas com genótipo "--" (não determinado)
- Database: verificar integridade, criar se não existir
- Genótipo: aceitar apenas ATCG (2 alelos)

**Output Esperado:**
- Tabela interativa com resultados
- Exemplo: 500 SNPs upload → encontra 120 matches em GWAS
- Pode ser exportado (no Nível 1, só visualização)

---

### NÍVEL 2: Impact Scoring + Busca Avançada

**Funcionalidade Adicional:**
1. Coluna "Impact Score" em cada resultado (1-10)
2. Ordenação padrão por Impact Score (descrescente = mais impactante)
3. Filtro slider: Impact Score mínimo
4. Filtro slider: P-value máximo (log scale)
5. Filtro dropdown: Categoria de trait (Metabolic, Cardiovascular, Neuropsych, Physical, Other)
6. Busca em texto livre: search em traits + genes
7. Todos os filtros aplicam em TEMPO REAL (sem botão submit)
8. Indicador visual: quantos resultados após filtros

**Características Técnicas Nível 2:**
- Full-text search (FTS5 SQLite) para performance
- Cálculo de score leva em conta p-value E allele frequency
- Fuzzy matching em busca (aceita "diabetes" para encontrar "Type 2 Diabetes")
- Cache de resultados para filtros repetidos

---

## ESPECIFICAÇÃO TÉCNICA DETALHADA

### Stack Tecnológico

```
LINGUAGEM: Python 3.9+
BACKEND:
  - sqlite3 (built-in, nenhuma dependência extra além das abaixo)
  - pandas (parsing 23andMe, data manipulation)
  - numpy (cálculos numéricos)

FRONTEND:
  - PyQt6 (desktop UI framework)
  
DATABASE:
  - SQLite (ficheiro gwas.db)
  - Pré-populado com GWAS Catalog + gnomAD allele frequencies
```

### Estrutura de Ficheiros Esperada

```
genetic_app/
├── main.py                              # Entrada aplicação
├── config.py                            # Constantes, paths, configurações
├── requirements.txt                     # Dependências pip
│
├── database/
│   ├── __init__.py
│   ├── setup_database.py                # Script criar + popular DB (com dados sample)
│   └── [gwas.db será gerado]
│
├── backend/
│   ├── __init__.py
│   ├── parsers.py                       # Parse ficheiro 23andMe
│   ├── search_engine.py                 # Query BD, matching logic, filtering
│   ├── scoring.py                       # Cálculo impact score
│   └── validators.py                    # Validação inputs
│
├── frontend/
│   ├── __init__.py
│   └── main_window.py                   # Interface principal (estrutura + lógica)
│
├── models/
│   ├── __init__.py
│   └── data_models.py                   # Dataclasses: SNPRecord, GWASMatch, FilterCriteria
│
├── utils/
│   ├── __init__.py
│   ├── logging_config.py                # Setup logging
│   └── file_utils.py                    # Utilitários ficheiros
│
└── tests/
    ├── __init__.py
    ├── test_parsers.py                  # Testes parsing
    ├── test_search_engine.py            # Testes matching + filtering
    └── test_scoring.py                  # Testes cálculo score
```

---

## DATABASE SCHEMA (SQLite)

### Tabela: gwas_variants

```
Coluna              | Tipo    | Descrição
────────────────────┼─────────┼─────────────────────────────────
id                  | INTEGER | Primary key
variant_id          | TEXT    | SNP identifier (rs6983267) [UNIQUE]
chromosome          | TEXT    | 1-22, X, Y, MT
position            | INTEGER | Genomic position
ref_allele          | TEXT    | Reference allele (A/T/C/G)
alt_allele          | TEXT    | Alternate allele
risk_allele         | TEXT    | Allele associated with higher risk
reported_trait      | TEXT    | Disease/trait name ("Colorectal cancer")
mapped_gene         | TEXT    | Gene symbol ("MYC")
p_value             | REAL    | GWAS p-value (5.2e-11)
odds_ratio          | REAL    | Effect size
sample_size         | INTEGER | Number of subjects in study
category            | TEXT    | Trait category
pubmed_id           | TEXT    | PubMed paper reference

Índices:
  - variant_id (PRIMARY)
  - p_value (for sorting)
  - category (for filtering)
  - reported_trait (for search)
```

### Tabela: allele_frequencies

```
Coluna              | Tipo    | Descrição
────────────────────┼─────────┼──────────────────────────────
id                  | INTEGER | Primary key
variant_id          | TEXT    | Foreign key to gwas_variants
af_overall          | REAL    | Allele frequency in all populations (0-1)
af_eur              | REAL    | European
af_afr              | REAL    | African
af_eas              | REAL    | East Asian
af_amr              | REAL    | American

Índices:
  - variant_id (FOREIGN KEY)
```

### Tabela Virtual: gwas_fts (Full-Text Search)

```
Tabela FTS5 para busca rápida em texto:
  - Indexa: variant_id, reported_trait, mapped_gene
  - Permite queries como: SELECT WHERE reported_trait MATCH "diabetes"
```

---

## FORMATO FICHEIRO 23andMe RAW DATA

### Exemplo de Input (genome_user.txt)

```
# Este é um ficheiro 23andMe raw data
# Format: RSID	CHROMOSOME	POSITION	GENOTYPE
#
# Linhas iniciadas com # são ignoradas
# RSID: SNP identifier começando com "rs"
# CHROMOSOME: 1-22, X, Y, MT
# POSITION: inteiro > 0
# GENOTYPE: 2 alelos (AA, AT, GG, etc.) ou "--" (não determinado)

rs3131972	1	694713	GG
rs12124819	1	713790	AG
rs11240777	1	856331	GG
rs6681049	1	909917	TT
rs4970383	1	1019440	AG
rs7899632	2	1234567	--
```

### Validação Parser

```
Para cada linha válida:
  ✓ Ignorar linhas vazias ou comentário
  ✓ RSID começa com "rs" + números
  ✓ CHROMOSOME em [1-22, X, Y, MT]
  ✓ POSITION é inteiro > 0
  ✓ GENOTYPE tem exatamente 2 caracteres de [ATCG--]
  ✓ Se genótipo é "--", skip linha
  ✓ Aceitar whitespace variable (tab, múltiplos espaços)
  ✗ Se erro: log warning com line number, continua processamento
  
Return: List[SNPRecord] ou exception se ficheiro invalido
```

---

## CÁLCULO DO IMPACT SCORE (Nível 2)

### Fórmula Completa

```
impact_score = score_p_value + score_allele_frequency

COMPONENTE 1: Score baseado em P-value
──────────────────────────────────────

score_p_value = 10 - min(10, log10(1 / p_value) / 10)

Explicação:
  - log10(1 / p_value) converte p-value em escala linear
  - Dividir por 10 normaliza (p_value=1e-10 → score ~9)
  - min(10, ...) caps em 10.0
  - Quanto MENOR p_value, MAIOR score

Exemplos:
  p_value = 1e-10:  log10(1e10) = 10.0, score = 10 - 1.0 = 9.0
  p_value = 1e-50:  log10(1e50) = 50.0, score = 10 - 5.0 = 5.0 (ou capped a 10)
  p_value = 1e-5:   log10(1e5)  = 5.0,  score = 10 - 0.5 = 9.5
  p_value = 0.001:  log10(1e3)  = 3.0,  score = 10 - 0.3 = 9.7
  p_value = 0.1:    log10(10)   = 1.0,  score = 10 - 0.1 = 9.9

COMPONENTE 2: Score baseado em Allele Frequency (raridade)
──────────────────────────────────────────────────────────

score_allele_frequency = (1 - af_overall) * 3

Explicação:
  - Variantes raras (AF baixa) são mais impactantes biologicamente
  - AF=0.5 (comum): score = 0.5 * 3 = 1.5
  - AF=0.1 (rara):  score = 0.9 * 3 = 2.7
  - AF=0.01 (muito rara): score = 0.99 * 3 = 2.97
  - Se AF não disponível: usar default AF=0.5

TOTAL:
──────

impact_score = min(10.0, max(0.0, score_p_value + score_allele_frequency))

Final clamp: resultado sempre em [0.0, 10.0]
```

### Implementação (Requisitos)

```
Função: calculate_impact_score(p_value: float, allele_frequency: float = None) -> float

Comportamento:
  - Input validation: p_value em (0, 1], AF em [0, 1] (ou None)
  - Edge cases:
    • p_value <= 0: raise ValueError
    • AF not in [0,1] and not None: raise ValueError
    • AF is None: use default 0.5
    • log10 error: log warning, return default score 5.0
  - Return: float in [0.0, 10.0]
  - Must be vectorizable (work with numpy arrays)
  
Example usage:
  score = calculate_impact_score(p_value=5.2e-11, allele_frequency=0.48)
  # Returns ~9.1

  scores = calculate_impact_score(
    p_value=np.array([1e-10, 1e-5, 0.01]),
    allele_frequency=0.5
  )
  # Returns array([9.0, 9.5, 9.7])
```

---

## OPERAÇÕES DE BUSCA E MATCHING

### 1. Matching User SNPs contra GWAS Catalog

```
Input:
  - List[SNPRecord] do ficheiro 23andMe (user_snps)
  - Conexão SQLite database (gwas_db)

Process:
  Para cada SNP user:
    1. Query: SELECT * FROM gwas_variants WHERE variant_id = ?
    2. Se encontrado (0+ rows):
       - Para cada row (possíveis múltiplas traits para mesmo SNP):
         a. Lookup allele frequency: SELECT af_overall FROM allele_frequencies
         b. Calcular impact_score usando formula acima
         c. Criar objeto GWASMatch com todos dados
    3. Se não encontrado: skip
  
Output:
  List[GWASMatch] com todos achados
  - Cada match contém: rsid, user_genotype, trait, gene, p_value, score, etc.
  - Um SNP user pode gerar múltiplos matches (múltiplos traits)

Performance requirement:
  - 500 SNPs user → query DB → ~100-200 matches em <2 segundos
```

### 2. Filtragem e Sorting (Nível 2)

```
Input:
  - List[GWASMatch] (todos resultados)
  - FilterCriteria object com:
    • min_score: float [0.0-10.0]
    • max_pvalue: float (e.g., 1e-5)
    • category: str ('Metabolic' | 'Cardiovascular' | 'ALL' | ...)
    • search_text: str (free text em trait + gene)
    • sort_by: 'score' | 'pvalue' | 'trait' (default: 'score')
    • sort_ascending: bool (default: False = descending)

Process:
  1. Filter score: manter apenas matches onde impact_score >= min_score
  2. Filter p-value: manter apenas matches onde p_value <= max_pvalue
  3. Filter category: se != 'ALL', manter apenas categoria selecionada
  4. Text search: match search_text em trait E gene (case-insensitive, fuzzy)
     - Usar FTS5 ou simple contains() com lower()
     - Exemplo: "diabetes" encontra "Type 2 Diabetes", "Diabetic neuropathy"
  5. Sort: ordenar por sort_key, ascending/descending conforme sort_ascending
  
Output:
  List[GWASMatch] filtrado e sorted

Real-time requirement:
  - Cada mudança filtro → re-apply filters + re-sort (sem delay visível)
  - Caching opcional para performance
```

### 3. Busca Full-Text (FTS5)

```
Configuração:
  CREATE VIRTUAL TABLE gwas_fts USING fts5(
    variant_id, reported_trait, mapped_gene,
    content=gwas_variants
  )

Usage:
  SELECT variant_id, reported_trait FROM gwas_fts 
  WHERE reported_trait MATCH 'diabetes'
  
Advantages:
  - Búsqueda rápida mesmo em 100k+ rows
  - Suporta boolean queries: "diabetes OR obesity"
  - Suporta prefix: "diab*"
  
Fallback (se performance aceitável):
  SELECT * FROM gwas_variants 
  WHERE reported_trait LIKE ? OR mapped_gene LIKE ?
```

---

## DATA MODELS (Dataclasses/Classes)

### SNPRecord

```
Representa um SNP lido do ficheiro 23andMe user.

Campos:
  - rsid: str (e.g., "rs6983267")
  - chromosome: str (e.g., "8")
  - position: int (e.g., 128413305)
  - genotype: str (e.g., "GG", "AA", "AT")

Validação no constructor:
  - rsid formato "rs" + digits
  - chromosome em [1-22, X, Y, MT]
  - position > 0
  - genotype exatamente 2 caracteres [ATCG]
  - Raise ValueError se inválido
```

### GWASMatch

```
Representa um match entre SNP user e entrada GWAS Catalog.

Campos:
  - rsid: str
  - chromosome: str
  - position: int
  - user_genotype: str (genótipo do utilizador)
  - gene: str | None (mapped_gene from GWAS)
  - trait: str (reported_trait)
  - risk_allele: str (qual alelo aumenta risco)
  - p_value: float
  - odds_ratio: float | None
  - sample_size: int | None
  - category: str (trait category)
  - allele_frequency: float (af_overall from gnomAD)
  - impact_score: float (calculated)

Métodos úteis:
  - __lt__, __gt__, etc. para sorting (enable por impact_score)
  - __repr__ para logging
```

### FilterCriteria

```
Armazena estado de filtros (Nível 2).

Campos:
  - min_score: float = 0.0
  - max_pvalue: float = 1.0
  - category: str = 'ALL'
  - search_text: str = ''
  - sort_by: str = 'score'  ('score' | 'pvalue' | 'trait')
  - sort_ascending: bool = False

Constructor:
  - Validar ranges (score [0-10], pvalue [0-1])
  - Validar sort_by em lista válida
```

---

## PROCESSAMENTO E THREADING

### Background Processing Requirement

```
Operação: Parse ficheiro 23andMe + Match contra GWAS + Calcular scores

Problema: Pode levar 2-5 segundos para 500+ SNPs
Solução: Usar threading para evitar UI freeze

Requisitos:
  - UI permanece responsiva durante processamento
  - Mostrar progress indicator (progress bar ou spinning indicator)
  - Permitir cancelar operação
  - Emitir signals ao completar (success, error, progress)

Implementação (strategy):
  - Criar worker thread (QThread em PyQt6)
  - Thread executa: parse → match → score
  - Thread emite pyqtSignal ao terminar ou erro
  - Main thread recebe signal e atualiza UI
```

---

## ERROR HANDLING E LOGGING

### Logging Strategy

```
Requirement:
  - Todo erro deve ser logged (file + console)
  - Debug info para troubleshooting
  - User-facing errors devem ser claros (não técnicos)

Setup:
  - Log file: logs/app.log (rotating, max 10MB, keep 5 backups)
  - Console: INFO level
  - File: DEBUG level
  - Format: timestamp | level | module | message

Exemplo log entry:
  2025-01-15 14:23:45,123 - backend.parsers - WARNING - Invalid genotype "--" at line 45
  2025-01-15 14:23:46,456 - frontend.main_window - INFO - Processing complete: 120 matches found
```

### Error Cases e Handling

```
PARSING ERRORS:
  - Invalid RSID format → skip line, log warning
  - Invalid chromosome → skip line, log warning
  - Missing field → skip line, log warning
  - Genotype "--" (unknown) → skip line, info level
  - File not readable → dialog error, return
  
DATABASE ERRORS:
  - Database file missing → auto-create via setup_database.py
  - Database corrupt → show dialog, suggest re-create
  - Query timeout → log warning, show timeout dialog
  - Connection error → show dialog, exit gracefully
  
CALCULATION ERRORS:
  - p_value <= 0 → log warning, use default score 5.0
  - AF outside [0,1] → log warning, use default AF 0.5
  - log10 error → log error, use fallback
  
UI ERRORS:
  - File dialog cancelled → do nothing (normal)
  - Table rendering 10k+ rows → use pagination (50/page)
  - Filter operation on empty results → show "No matches"
```

---

## PERFORMANCE REQUIREMENTS

### Query Performance

```
Expected performance:
  - Parse 500 SNPs: < 1 second
  - Match 500 SNPs vs 100k GWAS entries: < 2 seconds
  - Calculate 200 scores: < 0.5 seconds
  - Apply filters to 200 results: < 0.1 seconds (real-time)
  - Total end-to-end: < 5 seconds

Optimization strategies:
  - Index database (variant_id PRIMARY, p_value, category)
  - Parameterized SQL queries (prevent SQLi, improve execution plan)
  - Cache allele frequencies (frequent lookups)
  - Lazy evaluation of filters (only compute visible rows)
  - Use numpy for batch score calculations if possible
```

### Memory Management

```
Constraint: Single user machine (8GB+ RAM expected)

Optimization:
  - Load full database into memory on startup (gwas.db ~1GB)
  - Store results in memory during session
  - Cache filter operations (avoid re-computing)
  - Pagination: only render visible rows (not all 1000+)
```

---

## REQUIREMENTS.TXT

```
PyQt6==6.7.0
pandas==2.2.0
numpy==1.24.4
```

Installation:
```bash
pip install -r requirements.txt
```

---

## DATABASE SETUP

### Script: database/setup_database.py

Responsabilidades:
  1. Create tables (gwas_variants, allele_frequencies, gwas_fts)
  2. Create indexes on relevant columns
  3. Populate with sample GWAS data (200+ real variants from GWAS Catalog)
  4. Ensure data integrity (foreign keys, NOT NULL constraints)

Input data:
  - Use real GWAS Catalog entries (NHGRI-EBI GWAS Catalog dump)
  - Include diverse traits: Metabolic, Cardiovascular, Neuropsych, Physical
  - Includes SNPs with multiple trait associations
  - Include allele frequencies from gnomAD
  
Output:
  - File: database/gwas.db (SQLite3 format)
  - Size: ~1GB (for full dataset)
  - Can be distributed with app or generated on first run

Script debe ser idempotent:
  - Se database existe, opcionalmente drop e re-create
  - Or update existing data (append mode)
```

---

## TESTING REQUIREMENTS

### Unit Tests (tests/ directory)

```
test_parsers.py:
  - test_parse_valid_23andme_file()
  - test_parse_invalid_genotype_skipped()
  - test_parse_missing_fields_skipped()
  - test_parse_invalid_rsid_skipped()
  - test_parse_counts_correct()
  
test_scoring.py:
  - test_impact_score_valid_inputs()
  - test_impact_score_edge_cases(p_value=1e-100, af=0)
  - test_impact_score_clamped_to_10()
  - test_impact_score_with_missing_af()
  - test_impact_score_vectorized(numpy arrays)
  
test_search_engine.py:
  - test_match_user_snps_found()
  - test_match_user_snps_not_found()
  - test_filter_by_score()
  - test_filter_by_pvalue()
  - test_filter_by_category()
  - test_search_text_fuzzy()
  - test_sort_by_score_descending()
  - test_combined_filters()

Execution:
  python -m pytest tests/ -v
```

---

## DATA SAMPLE FOR TESTING

### Included Sample GWAS Data

```
Mínimo 200 variant entries incluindo:

Exemplo 1:
  variant_id: rs6983267
  chromosome: 8
  position: 128413305
  risk_allele: G
  reported_trait: Colorectal cancer
  mapped_gene: MYC
  p_value: 5.2e-11
  sample_size: 201518
  category: Oncology
  af_overall: 0.48

Exemplo 2:
  variant_id: rs1333049
  chromosome: 9
  position: 22125500
  risk_allele: C
  reported_trait: Cardiovascular disease
  mapped_gene: CDKN2A
  p_value: 4.8e-14
  sample_size: 256000
  category: Cardiovascular
  af_overall: 0.35

Exemplo 3 (múltiplas traits mesmo SNP):
  variant_id: rs7903146
  traits: Type 2 Diabetes, Fasting glucose, Insulin resistance
  p_values: 2.3e-36, 1.1e-20, 5.5e-15
  category: Metabolic

Data coverage:
  - 50+ traits Metabolic (diabetes, obesity, cholesterol)
  - 40+ traits Cardiovascular (MI, stroke, blood pressure)
  - 30+ traits Neuropsychiatric (depression, schizophrenia)
  - 40+ traits Physical (height, BMI, eye color)
  - 40+ other traits
```

---

## CONFIGURAÇÃO (config.py)

```
Constantes que devem estar em config.py:

DATABASE_PATH = 'database/gwas.db'
LOG_DIR = 'logs'
LOG_LEVEL = logging.DEBUG
LOG_FILE_SIZE = 10 * 1024 * 1024  # 10MB
LOG_BACKUP_COUNT = 5

RESULTS_PER_PAGE = 50
DEFAULT_ALLELE_FREQUENCY = 0.5
DEFAULT_IMPACT_SCORE = 5.0

VALID_CHROMOSOMES = ['1', '2', ..., '22', 'X', 'Y', 'MT']
TRAIT_CATEGORIES = ['Metabolic', 'Cardiovascular', 'Neuropsychiatric', 
                    'Physical Trait', 'Infectious', 'Immune', 'Other']

RSID_PATTERN = r'^rs\d+$'
GENOTYPE_PATTERN = r'^[ATCG]{2}$'
```

---

## SUMÁRIO DE FUNCIONALIDADES

| Feature | Nível 1 | Nível 2 | Status |
|---------|---------|---------|--------|
| Upload 23andMe | ✓ | ✓ | Required |
| Parse + Validate | ✓ | ✓ | Required |
| Match GWAS | ✓ | ✓ | Required |
| Table Display | ✓ | ✓ | Required |
| Sort by Column | ✓ | ✓ | Required |
| Pagination | ✓ | ✓ | Required |
| Impact Score Calc | ✓ | ✓ | Required |
| Score-based Sort | | ✓ | Required |
| Score Slider Filter | | ✓ | Required |
| P-value Filter | | ✓ | Required |
| Category Filter | | ✓ | Required |
| Text Search | | ✓ | Required |
| Real-time Filters | | ✓ | Required |

---

## INSTRUÇÕES FINAIS PARA O MODELO

### Liberdade de Design

```
O modelo DEVE fazer decisões autónomas sobre:

UI/Layout:
  - Como organizar componentes (layout top-to-bottom, side-by-side, panels, etc.)
  - Qual widget para cada função (slider vs text input para scores?)
  - Styling (cores, fonts, spacing - usar PyQt6 stylesheets)
  - How to present progress/results (tabela, gráfico, cards, etc.)
  
Architecture:
  - Padrões design (MVC, MVP, etc.)
  - Como estruturar signal/slots em PyQt6
  - Where to cache data (memory vs disk)
  - Threading model details
  
Code Organization:
  - Function naming (model pode escolher)
  - Class hierarchies (se preciso)
  - Module dependencies
  - Import organization

O modelo NÃO pode:
  - Criar pseudocode ou TODOs
  - Omitir implementação de qualquer feature
  - Usar dependências fora da lista
  - Violar os requisitos de performance
```

### Quality Expectations

```
✓ Complete, production-ready code (não tutorial/educational)
✓ All functions have docstrings (Args/Returns/Raises)
✓ All imports explicit (no 'import *')
✓ All error cases handled with try/except
✓ All user messages are clear (não técnico)
✓ All SQL parametrized (prevent injection)
✓ Type hints on all functions
✓ Logging at appropriate levels (DEBUG, INFO, WARNING, ERROR)
✓ Code is testable (dependencies injectable)
✓ No hardcoded paths (use config.py)
```

### Delivery Format

```
Generate:
  1. All files listed in "Estrutura de Ficheiros"
  2. requirements.txt with exact versions
  3. database/setup_database.py with sample data (ready to run)
  4. Complete working code (no stubs)
  5. tests/ directory with full coverage
  6. README or docstring explaining how to run

Output should be:
  - Organized by file (markdown code blocks with filename)
  - Ready to copy-paste into respective files
  - Runnable as-is after: pip install -r requirements.txt && python main.py
```

---

## END OF SPECIFICATION

Este prompt fornece:
- ✓ Especificação funcional completa (Níveis 1-2)
- ✓ Requisitos técnicos detalhados (BD, algoritmos, performance)
- ✓ Arquitetura abstrata (sem pseudocode)
- ✓ Liberdade de design UI/código ao modelo
- ✓ Criterios de qualidade e testes
- ✓ Dados sample e scripts setup

Pronto para ser enviado diretamente a Claude Opus 4.5 ou modelo equivalente.
