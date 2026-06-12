#!/usr/bin/env python3

import os
import re
import gzip
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox
import matplotlib.pyplot as plt
from collections import Counter

# Sequenze per analisi interruzione CAA (LOSS globale vicino a CCGCCA)
WT_SEQ = "CAACAGCCGCCA"      # con interruzione CAA
LOSS_SEQ = "CAGCAGCCGCCA"    # perdita dell'interruzione CAA

# Sequenze per analisi GAIN della CAA in 8 posizioni (come da specifiche)
GAIN_POSITIONS = [
    (1, "CAGCAACAGCCGCCA",                     "CAACAACAGCCGCCA"),
    (2, "CAGCAGCAACAGCCGCCA",                  "CAACAGCAACAGCCGCCA"),
    (3, "CAGCAGCAGCAACAGCCGCCA",               "CAACAGCAGCAACAGCCGCCA"),
    (4, "CAGCAGCAGCAGCAACAGCCGCCA",            "CAACAGCAGCAGCAACAGCCGCCA"),
    (5, "CAGCAGCAGCAGCAGCAACAGCCGCCA",         "CAACAGCAGCAGCAGCAACAGCCGCCA"),
    (6, "CAGCAGCAGCAGCAGCAGCAACAGCCGCCA",      "CAACAGCAGCAGCAGCAGCAACAGCCGCCA"),
    (7, "CAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCA",   "CAACAGCAGCAGCAGCAGCAGCAACAGCCGCCA"),
    (8, "CAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCA","CAACAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCA"),
]

def count_repeats_fastq(fastq_paths, motif, min_len, min_repeats=5):
    """
    - Conta le ripetizioni STR (CAG) con filtri e genera istogramma (frequenze >=25).
    - LOSS: valuta WT vs LOSS globali, scrive TXT con WT/(WT+LOSS)*100 e interpretazione;
            se 'assente', aggiunge distribuzione CAG nelle reads LOSS.
    - GAIN: per ciascuna delle 8 posizioni, valuta WT_pos vs GAIN_pos, calcola WT/(WT+GAIN)*100,
            scrive interpretazione; se 'presente', aggiunge distribuzione CAG nelle reads GAIN_pos.
    - Salva PNG automaticamente accanto al TSV:
        - basename.png (istogramma principale)
        - *_CAA_summary.png (unico plot % per LOSS + GAIN1..8, etichette sulla barra con % e N)
        - *_CAA_distributions.png (9 pannelli: distribuzioni CAG per LOSS e GAIN1..8)
    - Un solo popup finale elenca tutti i file salvati.
    """
    motif_esc = re.escape(motif)
    motif_len = len(motif)
    pattern = re.compile(f'(?:{motif_esc})+')

    counts = []

    # Variabili LOSS globale
    total_reads = 0
    wt_reads = 0
    loss_reads = 0
    loss_repeat_counts = []

    # Variabili GAIN per posizioni 1..8
    wtpos_counts = {pos: 0 for (pos, _, _) in GAIN_POSITIONS}
    gain_counts = {pos: 0 for (pos, _, _) in GAIN_POSITIONS}
    gain_repeat_counts = {pos: [] for (pos, _, _) in GAIN_POSITIONS}

    # Per riepilogo salvataggi (un solo OK finale)
    saved_files = []

    # Lettura FASTQ / FASTQ.GZ
    for fastq_path in fastq_paths:
        try:
            fh = gzip.open(fastq_path, 'rt') if fastq_path.lower().endswith('.gz') else open(fastq_path, 'r')
        except Exception as e:
            messagebox.showerror('Error', f'Cannot open FASTQ file {fastq_path}: {e}')
            continue

        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline(); fh.readline()  # skip '+' e qualità
            total_reads += 1

            seq_u = seq.upper()

            # Calcolo una volta il max run contigua di CAG sulla read
            runs_all = pattern.findall(seq_u)
            abs_count_all = (len(max(runs_all, key=len)) // motif_len) if runs_all else 0

            # --- LOSS globale ---
            has_wt = WT_SEQ in seq_u
            has_loss = LOSS_SEQ in seq_u
            if has_wt:
                wt_reads += 1
            if has_loss:
                loss_reads += 1
                loss_repeat_counts.append(abs_count_all)

            # --- GAIN posizionali (1..8) ---
            for (pos, wtpos, gainpos) in GAIN_POSITIONS:
                if wtpos in seq_u:
                    wtpos_counts[pos] += 1
                if gainpos in seq_u:
                    gain_counts[pos] += 1
                    gain_repeat_counts[pos].append(abs_count_all)

            # --- Parte "conta triplette" con filtri per TSV/istogramma principale ---
            if len(seq) < min_len:
                continue
            runs = pattern.findall(seq_u)
            abs_count = (len(max(runs, key=len)) // motif_len) if runs else 0
            if abs_count < min_repeats:
                continue
            counts.append(abs_count)

        fh.close()

    # Se nessuna read supera i filtri per la parte "triplette", come l'originale
    if not counts:
        messagebox.showwarning('No Data', 'No reads passed filters; nothing to plot.')
        return

    # Salva TSV dei conteggi
    save_tsv = filedialog.asksaveasfilename(
        title='Save repeat counts TSV',
        defaultextension='.tsv',
        filetypes=[('TSV files','*.tsv')]
    )
    if not save_tsv:
        messagebox.showwarning('Cancelled','No TSV file specified. Outputs not saved.')
        return

    with open(save_tsv, 'w') as out:
        out.write('read_name\tabs_repeats\n')
        for i, c in enumerate(counts, start=1):
            out.write(f'read_{i}\t{c}\n')
    saved_files.append(f"TSV\t{save_tsv}")

    # Percorsi di output
    base_noext, _ = os.path.splitext(save_tsv)
    png_main_path   = f"{base_noext}.png"
    txt_path        = f"{base_noext}_CAA_interruption.txt"
    png_summary     = f"{base_noext}_CAA_summary.png"
    png_distribs    = f"{base_noext}_CAA_distributions.png"

    # --- TXT: LOSS globale ---
    denom_loss = wt_reads + loss_reads
    if denom_loss == 0:
        ratio_loss = None
        status_loss = "Interruzione CAA non valutabile (nessuna read WT o LOSS trovata)"
    else:
        ratio_loss = (wt_reads / denom_loss) * 100.0
        if ratio_loss > 90:
            status_loss = "Interruzione CAA presente"
        elif ratio_loss < 90:
            status_loss = "Interruzione CAA assente"
        else:
            status_loss = "Interruzione CAA borderline (rapporto = 90)"

    with open(txt_path, 'w') as fout:
        fout.write("# Analisi interruzione CAA vicino alle triplette CAG (LOSS globale)\n")
        fout.write(f"# FASTQ analizzati: {len(fastq_paths)}\n")
        fout.write(f"Total_reads\t{total_reads}\n")
        fout.write(f"WT_sequence\t{WT_SEQ}\n")
        fout.write(f"WT_count_reads\t{wt_reads}\n")
        fout.write(f"LOSS_sequence\t{LOSS_SEQ}\n")
        fout.write(f"LOSS_count_reads\t{loss_reads}\n")
        fout.write("WT_to_LOSS_ratio_percent\t" + ("NA" if ratio_loss is None else f"{ratio_loss:.3f}") + "\n")
        fout.write(f"Interpretation\t{status_loss}\n")

        if denom_loss > 0 and ratio_loss is not None and ratio_loss < 90 and loss_reads > 0:
            fout.write("\n# Distribuzione ripetizioni CAG nelle reads con assenza di CAA (LOSS)\n")
            dist = Counter(loss_repeat_counts)
            for rep in sorted(dist.keys()):
                fout.write(f"CAG_repeats={rep}\treads={dist[rep]}\n")

        # --- TXT: GAIN posizionali 1..8 ---
        fout.write("\n# Analisi GAIN di CAA in posizioni specifiche (1..8)\n")
        for (pos, wtpos, gainpos) in GAIN_POSITIONS:
            wt_c = wtpos_counts[pos]
            g_c = gain_counts[pos]
            denom_pos = wt_c + g_c
            if denom_pos == 0:
                ratio_pos = None
                status_pos = f"Guadagno CAA posizione {pos} non valutabile (nessuna read WT_pos o GAIN_pos trovata)"
            else:
                ratio_pos = (wt_c / denom_pos) * 100.0  # WT/(WT+GAIN)*100
                if ratio_pos > 90:
                    status_pos = f"Guadagno CAA posizione {pos} assente"
                elif ratio_pos < 90:
                    status_pos = f"Guadagno CAA posizione {pos} presente"
                else:
                    status_pos = f"Guadagno CAA posizione {pos} borderline (rapporto = 90)"

            fout.write(f"\nPOS\t{pos}\n")
            fout.write(f"WT_pos_sequence\t{wtpos}\n")
            fout.write(f"GAIN_pos_sequence\t{gainpos}\n")
            fout.write(f"WT_pos_count_reads\t{wt_c}\n")
            fout.write(f"GAIN_pos_count_reads\t{g_c}\n")
            fout.write("WTpos_to_GAINpos_ratio_percent\t" + ("NA" if ratio_pos is None else f"{ratio_pos:.3f}") + "\n")
            fout.write(f"Interpretation\t{status_pos}\n")

            # Se guadagno presente, distribuzione CAG nelle reads GAIN_pos
            if denom_pos > 0 and ratio_pos is not None and ratio_pos < 90 and g_c > 0:
                fout.write("# Distribuzione CAG (max run) nelle reads GAIN_pos\n")
                dist_g = Counter(gain_repeat_counts[pos])
                for rep in sorted(dist_g.keys()):
                    fout.write(f"CAG_repeats={rep}\treads={dist_g[rep]}\n")

    saved_files.append(f"TXT\t{txt_path}")

    # --- Istogramma principale (conte ripetizioni filtrate, come prima) ---
    freq = Counter(counts)
    filtered = {k: v for k, v in freq.items() if v >= 25}
    if filtered:
        x = sorted(filtered.keys())
        y = [filtered[k] for k in x]
        fig, ax = plt.subplots(figsize=(10, 6))
        min_x, max_x = x[0], x[-1]
        ax.axvspan(37.5, 39.5, facecolor='red', alpha=0.1)
        ax.axvspan(39.5, max_x + 0.5, facecolor='red', alpha=0.2)
        bars = ax.bar(x, y, color='blue')
        ax.set_xlabel('Number of repeats')
        ax.set_ylabel('Read count')
        ax.set_title(f'Histogram of {motif} repeats (counts ≥25)')
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + max(y)*0.01,
                    f'{int(height)}', ha='center', va='bottom', fontsize=7, color='gray')
        ax.set_xticks(x)
        ax.set_xticklabels([str(val) for val in x], fontsize=8)
        plt.tight_layout()
        fig.savefig(png_main_path)
        plt.close(fig)
        saved_files.append(f"PNG (hist)\t{png_main_path}")
    else:
        # nessun istogramma principale
        pass

    # =========================
    #  PLOT UNICO: LOSS + GAIN
    # =========================
    positions = [pos for (pos, _, _) in GAIN_POSITIONS]

    # percentuali e denominatori
    def ratio_and_denom(w, x):
        denom = w + x
        if denom == 0:
            return None, 0
        return (w / denom) * 100.0, denom

    ratio_loss_val, denom_loss_val = ratio_and_denom(wt_reads, loss_reads)
    ratio_gain_vals = []
    denom_gain_vals = []
    for pos in positions:
        r, d = ratio_and_denom(wtpos_counts[pos], gain_counts[pos])
        ratio_gain_vals.append(r)
        denom_gain_vals.append(d)

    cats = ["LOSS"] + [f"G{p}" for p in positions]
    ratios = [ratio_loss_val] + ratio_gain_vals
    denoms = [denom_loss_val] + denom_gain_vals
    ys = [0 if r is None else r for r in ratios]

    figS, axS = plt.subplots(figsize=(12, 6))
    bars = axS.bar(range(len(cats)), ys)
    axS.axhline(90, linestyle='--', linewidth=1)
    axS.set_xticks(range(len(cats)))
    axS.set_xticklabels(cats)
    axS.set_ylim(0, 100)
    axS.set_ylabel('WT / (WT + X) %')
    # titolo un po' più in alto
    axS.set_title('CAA: LOSS & GAIN summary', pad=24)

    # Etichette SULLA BARRA: "%\n(N=reads_usate)"
    for i, (bar, r, n) in enumerate(zip(bars, ratios, denoms)):
        if r is None:
            label = "NA\n(N=0)"
            y_pos = bar.get_height() + 2  # appena sopra la barra vuota
            axS.text(bar.get_x() + bar.get_width()/2, y_pos, label,
                     ha='center', va='bottom', fontsize=9)
        else:
            label = f"{r:.1f}%\n(N={n})"
            h = bar.get_height()
            # se la barra è piccola, scrivo appena sopra; altrimenti dentro
            if h < 8:
                y_pos = h + 2
                axS.text(bar.get_x() + bar.get_width()/2, y_pos, label,
                         ha='center', va='bottom', fontsize=9)
            else:
                axS.text(bar.get_x() + bar.get_width()/2, h * 0.5, label,
                         ha='center', va='center', fontsize=9)

    plt.tight_layout()
    figS.savefig(png_summary)
    plt.close(figS)
    saved_files.append(f"PNG (summary)\t{png_summary}")

    # ===============================================
    #  FILE UNICO: TUTTE LE DISTRIBUZIONI (LOSS+GAIN)
    # ===============================================
    all_titles = [f"LOSS (N={loss_reads})"] + [f"G{p} (N={gain_counts[p]})" for p in positions]
    all_distribs = [Counter(loss_repeat_counts)] + [Counter(gain_repeat_counts[p]) for p in positions]

    figD, axes = plt.subplots(3, 3, figsize=(14, 12))
    axes = axes.flatten()

    for idx, (title, dist) in enumerate(zip(all_titles, all_distribs)):
        axd = axes[idx]
        if dist:
            x_vals = sorted(dist.keys())
            y_vals = [dist[k] for k in x_vals]
            axd.bar(x_vals, y_vals)
            axd.axvspan(37.5, 39.5, facecolor='red', alpha=0.1)
            axd.axvspan(39.5, (max(x_vals) + 0.5), facecolor='red', alpha=0.2)
            axd.set_xlabel('CAG repeats')
            axd.set_ylabel('Read count')
            axd.set_title(title, fontsize=10)
        else:
            axd.axis('off')
            axd.set_title(title, fontsize=10)
            axd.text(0.5, 0.5, 'Nessuna read', ha='center', va='center', fontsize=10)
    plt.tight_layout()
    figD.savefig(png_distribs)
    plt.close(figD)
    saved_files.append(f"PNG (distributions)\t{png_distribs}")

    # Unico popup finale con l’elenco dei file salvati
    if saved_files:
        msg = "Sono stati generati i seguenti file:\n" + "\n".join(f"- {p}" for p in saved_files)
        messagebox.showinfo('Salvataggio completato', msg)


if __name__ == '__main__':
    root = tk.Tk(); root.withdraw()
    messagebox.showinfo('Input Required', 'Select FASTQ/FASTQ.GZ files and parameters')

    fastq_paths = filedialog.askopenfilenames(
        title='Select FASTQ or FASTQ.GZ files',
        filetypes=[('FASTQ/FASTQ.GZ','*.fastq *.fq *.fastq.gz *.fq.gz')]
    )
    if not fastq_paths:
        messagebox.showwarning('Cancelled','No FASTQ files selected. Exiting.'); exit()

    motif = simpledialog.askstring('Motif','Enter repeat motif (e.g. CAG):',initialvalue='CAG')
    if not motif:
        messagebox.showwarning('Cancelled','No motif provided. Exiting.'); exit()

    min_len = simpledialog.askinteger('Min Read Length','Minimum read length to include:',initialvalue=500, minvalue=1)
    count_repeats_fastq(list(fastq_paths), motif, min_len)
