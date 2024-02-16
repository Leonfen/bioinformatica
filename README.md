# Progetto di Bioinformatica

Questo progetto prevede l'analisi dei dati sull'espressione genica utilizzando vari test statistici e tecniche di visualizzazione. Di seguito sono riportate le principali funzioni implementate nello script `main.py`.

## Panoramica delle Funzioni

### Recupero e Pre-elaborazione dei Dati

1. **`main()`**
   - Punto di ingresso dello script.
   - Recupera i dati sull'espressione genica utilizzando il modulo GEO.
   - Esegue passaggi di pre-elaborazione come l'assegnazione di classi e la suddivisione dei dati in gruppi di controllo e infetti.

### Analisi Statistica

2. **`t_test()`**
   - Esegue un test t per confrontare l'espressione genica tra gruppi di controllo e infetti.
   - Restituisce valori p e nomi dei geni.

3. **`oneway()`**
   - Esegue un test ANOVA unidirezionale per confrontare l'espressione genica tra gruppi.
   - Restituisce valori p e nomi dei geni.

4. **`kruskal()`**
   - Esegue un test di Kruskal-Wallis H per confrontare l'espressione genica tra gruppi.
   - Restituisce valori p e nomi dei geni.

5. **`utest()`**
   - Esegue un test U di Mann-Whitney per confrontare l'espressione genica tra gruppi.
   - Restituisce valori p e nomi dei geni.

### Visualizzazione

6. **`plots.plotting()`**
   - Genera grafici Volcano per visualizzare le variazioni nell'espressione genica.
   - Accetta risultati statistici e dataframe.
   - Plotta i geni differenzialmente espressi.

### Output

7. **`utils.output.all_gene_names()`**
   - Restituisce i nomi dei geni categorizzati come sovraespressi e sottopressi.

## Utilizzo

Per eseguire l'analisi e la visualizzazione, eseguire la funzione `main()`.

```python
python main.py
```

## Dipendenze

- `scipy`
- Moduli personalizzati: `utils`, `graph`

Assicurarsi di avere installate le dipendenze richieste per eseguire lo script correttamente.

## Nota

Regolare le soglie di significatività e i parametri in base al design sperimentale e ai requisiti specifici.