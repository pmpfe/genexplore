# PROMPT REVISTO: GenExplore - Integração de Análise Poligénica

## OBJETIVO GLOBAL

Refatorar a aplicação GenExplore para evoluir de uma ferramenta de análise de SNPs individuais para uma plataforma de análise genética compreensiva.

**Scope:**
1. Manter 100% das features monogénicas existentes (sem alterações funcionais)
2. Adicionar novo módulo para análise de risco poligénico
3. Implementar interface unificada usando tabs
4. Adicionar capacidade de controlo de versões de dados

---

## ESCOPO FUNCIONAL DETALHADO

### TAB 1: ANÁLISE MONOGÉNICA (Existente - Refatorado)

**O que está:** Funcionalidade atual de GenExplore
**O que muda:** Apenas a integração na UI (move para primeiro tab)
**Preservação:** 100% das features e comportamentos existentes

- Upload de ficheiro genético (23andMe raw data)
- Busca e matching em base de dados GWAS
- Cálculo de scores de impacto individual
- Tabela interativa com filtragem e ordenação
- Paginação de resultados

---

### TAB 2: ANÁLISE POLIGÉNICA (Novo)

#### Feature 1: Navegador de Scores Poligénicos

**Vista Principal:** Tabela de scores disponíveis

Elementos visíveis ao utilizador:
- Nome da trait (doença ou característica genética)
- Categoria da trait (ex: Metabolismo, Cardiovascular, etc.)
- Informação sobre o estudo (ano de publicação, população estudada)
- Número de variantes envolvidas
- [Se pré-computado] Score e percentil do utilizador
- [Se pré-computado] Categoria de risco (Baixo/Intermédio/Alto)
- Botão para análise detalhada

**Interação:**
- Filtros por categoria
- Busca textual em nomes de traits
- Ordenação por qualquer coluna
- Sem carregamento por página (todos scores sempre disponíveis, mesmo que 5000+)

---

#### Feature 2: Análise Detalhada de Um Score

**O que acontece ao utilizador clicar num score:**

1. **Computação:**
   - Calcular risco poligénico do utilizador para esse score específico
   - Usar genótipo do ficheiro carregado (monogenic tab)
   - Lidar graciosamente com SNPs ausentes

2. **Visualização do Resultado:**
   - Mostrar distribuição da população de referência
   - Marcar claramente onde o utilizador se situa nessa distribuição
   - Indicar percentil (ex: "85º percentil")
   - Indicar categoria de risco (Baixo/Intermédio/Alto)
   - Mostrar score bruto e normalizado

3. **Contexto Científico:**
   - Detalhes do estudo (autores, publicação, DOI)
   - População de referência usada
   - Número de participantes no estudo
   - Disclaimer adequado ("Ferramenta de pesquisa, não diagnóstico médico")

4. **Informação de Qualidade:**
   - Percentagem de variantes encontradas no genótipo do utilizador
   - Aviso se cobertura for baixa
   - Contexto sobre aplicabilidade a população do utilizador (se diferente da referência)

5. **Exploração:**
   - Opção de ver quais variantes individuais contribuem mais para o score
   - Opção de exportar resultados (formato simples, ex: PDF ou texto)

---

#### Feature 3: Gestão de Bases de Dados

**Interface de Configuração/Definições**

Informação Visível:
- Versão atual da base de dados GWAS Catalog
- Data de última atualização
- Versão atual da base de dados de Scores Poligénicos
- Data de última atualização

Capacidades:
- Botão para "Verificar atualizações disponíveis"
- Botão para "Atualizar agora" (se atualização disponível)
- Indicador de progresso durante atualização
- Confirmação após conclusão bem-sucedida

Comportamento de Atualização:
- Descarregar novo ficheiro de base de dados completo de fonte oficial
- Validar integridade do ficheiro descarregado
- Manter cópia de segurança da versão anterior
- Ativar nova base de dados atomicamente (tudo-ou-nada, sem intermédio)
- Registar update em log (para auditoria)

---

### TAB 0: Componentes Partilhados

**Upload de Genótipo:**
- Interface unificada para carregar ficheiro genético
- Usada por ambos tabs (monogénico e poligénico)
- Upload acontece uma vez, resultados disponíveis em ambos tabs

**Barra de Status:**
- Indica se ficheiro genético está carregado
- Mostra informação básica sobre genótipo (número de SNPs, cobertura geral)

---

## CONCEITOS-CHAVE SEM ESPECIFICAÇÃO TÉCNICA

### Risco Poligénico

**O que é:**
- Combinação ponderada de efeitos genéticos de múltiplas variantes
- Cada variante tem um "peso" (quanto contribui para o risco)
- Soma de contribuições = risco do indivíduo

**Interpretação:**
- Convertido para escala poblacional para ser interpretável
- Posição do utilizador relativa à população (percentil)
- Categorização em níveis (baixo, intermédio, alto)

**Limitações a Comunicar:**
- Nem todos os SNPs estarão presentes no genótipo
- Importância de reportar cobertura atingida
- Variabilidade por ancestry (score pode ser mais/menos acurado dependendo de população)
- Não é diagnóstico, apenas indicador de risco relativo

---

### Pré-computação de Scores

**Conceito:**
- Ao utilizador carregar genótipo, computar risco para TODOS os scores disponíveis
- Armazenar resultados em memória durante sessão
- Permitir exploração rápida (sem esperas adicionais)

**Timing:**
- Executar em background enquanto utilizador explora outras features
- Mostrar progresso
- Atualizar tabela de scores conforme resultados ficam disponíveis

**Benefício:**
- Navegação instantânea (sem recomputação ao clicar)
- Permite ordenação/filtragem por scores do utilizador
- Melhor experiência de utilizador

---

### Atualização de Bases de Dados

**Modelo Offline + Batch:**
- Todas as bases de dados residem localmente
- Atualização feita em batch (download completo, não incremental)
- Processo atómico (tudo funciona ou nada muda)
- Sem dependência de internet durante uso (apenas na atualização)

**Vantagens:**
- Simplicidade de implementação
- Integridade garantida (sem múltiplas versões parciais)
- Fácil de fazer backup
- Fácil de rollback se algo correr mal

---

## REQUISITOS DE QUALIDADE

### Funcionalidade
- ✅ Todos os cálculos devem poder ser validados contra ferramentas de referência
- ✅ Lidar graciosamente com dados incompletos (SNPs ausentes)
- ✅ Mensagens de erro em linguagem de utilizador
- ✅ Documentação clara de limitações e assumptions

### Performance
- Cálculo de um score: <100 ms
- Pré-computação completa: <2 minutos (aceitável em background)
- UI permanece responsiva durante computação
- Carregamento de visualizações: <200 ms

### Segurança e Confiabilidade
- Dados genéticos nunca deixam o computador
- Backup automático antes de atualizar BD
- Validação de ficheiros descarregados
- Recuperação graceful de falhas
- Nenhuma perda de dados do utilizador

### Logging e Debugging
- Logs claros para troubleshooting
- Rastreamento de updates de BD
- Informação de erros útil para debugging
- Sem informação sensível em logs públicos

---

## ARQUITETURA GERAL (Apenas Conceitual)

### Camada Apresentação
- Interface tabbed unificada
- Tab 1: Análise monogénica
- Tab 2: Análise poligénica com navegador e visualizador
- Componentes partilhados (upload, status)

### Camada Lógica
- Motor de cálculo monogénico (existente)
- Motor de cálculo poligénico (novo)
- Gestor de atualização de bases de dados
- Sistema de caching de resultados pré-computados
- Handlers de erro consistentes

### Camada Dados
- Base de dados monogénica (GWAS Catalog)
- Base de dados poligénica (Scores e metadados)
- Tabela de versioning para rastreamento de atualizações
- Capacidade de suportar múltiplos ancestry/populações

---

## LIBERDADE DE DESIGN

O desenvolvedor tem total liberdade para:

**Implementação:**
- Escolher estrutura de código e organização
- Decidir padrões de design (MVC, etc.)
- Selecionar bibliotecas para visualização
- Definir estratégia de threading/background processing
- Determinar formato interno de caching

**Integração:**
- Como estruturar comunicação entre tabs
- Como gerenciar estado compartilhado
- Como atualizar UI progressivamente

**User Experience:**
- Design específico de interfaces
- Layout de componentes
- Escolha de cores e tipografia
- Detalhes de animações ou transições

---

## RESTRIÇÕES FIXAS

**O que NÃO pode mudar:**
- ✅ Funcionalidade monogénica mantém 100% compatibilidade
- ✅ Offline-first (sem requisito de internet para uso)
- ✅ Atualização de BD em batch completo
- ✅ Pré-computação possível sem latência crítica
- ✅ Cálculos validáveis contra referência
- ✅ Privacidade (genotypes local-only)

---

## REQUISITOS DE ENTREGA

**Estrutura do Código:**
- Modular e bem organizado
- Completamente funcional (sem TODOs ou placeholders)
- Type hints em funções públicas
- Docstrings em classes e funções principais
- Error handling robusto

**Testes:**
- Testes unitários para lógica de cálculo crítica
- Testes de integração para fluxos principais
- Testes de performance (verificar requisitos de timing)
- Dados sample para validação

**Documentação:**
- Como construir e executar aplicação
- Como usar a interface
- Limitações conhecidas
- Como validar resultados contra ferramentas de referência

**Qualidade:**
- Código limpo e legível
- Sem hardcoding de valores
- Logging apropriado em múltiplos níveis
- Tratamento graceful de erros

---

## PRÓXIMOS PASSOS APÓS IMPLEMENTAÇÃO

1. Validação contra ferramentas oficiais de scoring
2. Testing com dados reais de múltiplas populações
3. Otimizações de performance se necessário
4. Interface polishing baseado em feedback de utilizadores
5. Documentação para utilizadores finais

---

## NOTAS IMPORTANTES

**Sobre Dados Genéticos:**
- Cobertura: nem todos SNPs estarão em cada score
- Ancestry: scores são calibrados para populações específicas, accuracy varia
- Assembly: genoma pode estar em versões diferentes, pode afetar matching

**Sobre Visualização:**
- Distribuição populacional é forma standard de comunicar risco
- Percentil é mais intuitivo que score bruto para utilizadores leigos
- Cores devem ser accessível para utilizadores com daltonismo

**Sobre Performance:**
- Pré-computação em background é preferível a esperas na UI
- Caching em memória aceita porque sessão é típicamente <1 hora

---

## STATUS

✅ Especificação Funcional Completa
✅ Requisitos de Qualidade Definidos
✅ Liberdade de Design Preservada
✅ Pronto para Implementação

