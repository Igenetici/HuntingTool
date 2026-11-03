"""Graphical interface (Tkinter, included with standard Python installations)."""
from __future__ import annotations

import os
import queue
import threading
import traceback
import webbrowser
from pathlib import Path

from . import __version__


def main(inputs: list[str] | None = None, outdir: str | None = None,
         autorun: bool = False, autoclose_ms: int = 0) -> int:
    """Open the GUI. The optional arguments pre-fill the form (used by tests and launchers)."""
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except ImportError:
        print("Tkinter is not available. Install it (e.g. 'sudo apt install python3-tk') "
              "or use the command line: huntingtool run -i INPUT -o OUTDIR")
        return 1
    from .fastq import discover_samples
    from .pipeline import Options, run

    root = tk.Tk()
    root.title(f"HuntingTool {__version__}")
    root.geometry("860x640")
    root.minsize(720, 520)
    q: queue.Queue = queue.Queue()
    state = {"report": None, "running": False}

    frm = ttk.Frame(root, padding=12)
    frm.pack(fill="both", expand=True)
    frm.columnconfigure(0, weight=1)

    # inputs
    ttk.Label(frm, text="Input FASTQ files or folders (a MinKNOW fastq_pass folder gives one sample per barcode)").grid(
        row=0, column=0, sticky="w")
    lst = tk.Listbox(frm, height=7, selectmode="extended")
    lst.grid(row=1, column=0, sticky="nsew", pady=(2, 4))
    bt = ttk.Frame(frm)
    bt.grid(row=1, column=1, sticky="n", padx=(8, 0))

    def add_files():
        for f in filedialog.askopenfilenames(title="Select FASTQ files",
                                             filetypes=[("FASTQ", "*.fastq *.fq *.fastq.gz *.fq.gz"), ("All", "*")]):
            lst.insert("end", f)

    def add_folder():
        d = filedialog.askdirectory(title="Select a folder")
        if d:
            lst.insert("end", d)

    def remove():
        for i in reversed(lst.curselection()):
            lst.delete(i)

    for txt, cmd in (("Add files…", add_files), ("Add folder…", add_folder), ("Remove", remove),
                     ("Clear", lambda: lst.delete(0, "end"))):
        ttk.Button(bt, text=txt, command=cmd, width=14).pack(fill="x", pady=2)

    # output
    ttk.Label(frm, text="Output folder").grid(row=2, column=0, sticky="w", pady=(6, 0))
    outv = tk.StringVar(value=str(Path.home() / "HuntingTool_results"))
    of = ttk.Frame(frm)
    of.grid(row=3, column=0, columnspan=2, sticky="ew")
    of.columnconfigure(0, weight=1)
    ttk.Entry(of, textvariable=outv).grid(row=0, column=0, sticky="ew")
    ttk.Button(of, text="Browse…", command=lambda: outv.set(filedialog.askdirectory() or outv.get())).grid(
        row=0, column=1, padx=(6, 0))

    # parameters
    pf = ttk.LabelFrame(frm, text="Parameters", padding=8)
    pf.grid(row=4, column=0, columnspan=2, sticky="ew", pady=8)
    modev = tk.StringVar(value="htt")
    motifv = tk.StringVar(value="CAG")
    minlen = tk.IntVar(value=100)
    minq = tk.DoubleVar(value=0.0)
    minfrac = tk.DoubleVar(value=0.10)
    thr = tk.IntVar(value=max(1, (os.cpu_count() or 2) - 1))
    wreads = tk.BooleanVar(value=False)
    merge = tk.BooleanVar(value=False)
    fields = [("Mode", ttk.Combobox(pf, textvariable=modev, values=["htt", "motif"], width=8, state="readonly")),
              ("Motif (motif mode)", ttk.Entry(pf, textvariable=motifv, width=8)),
              ("Min read length (bp)", ttk.Spinbox(pf, from_=0, to=100000, textvariable=minlen, width=8)),
              ("Min mean read Q", ttk.Spinbox(pf, from_=0, to=60, increment=1, textvariable=minq, width=8)),
              ("Min allele fraction", ttk.Spinbox(pf, from_=0.01, to=0.5, increment=0.01, textvariable=minfrac, width=8)),
              ("Parallel samples", ttk.Spinbox(pf, from_=1, to=64, textvariable=thr, width=8))]
    for i, (lab, w) in enumerate(fields):
        ttk.Label(pf, text=lab).grid(row=i // 3, column=(i % 3) * 2, sticky="w", padx=(0, 4), pady=2)
        w.grid(row=i // 3, column=(i % 3) * 2 + 1, sticky="w", padx=(0, 18), pady=2)
    ttk.Checkbutton(pf, text="Write per-read table", variable=wreads).grid(row=2, column=0, columnspan=2, sticky="w")
    ttk.Checkbutton(pf, text="Merge all inputs into one sample", variable=merge).grid(row=2, column=2, columnspan=3, sticky="w")

    # run + log
    rb = ttk.Frame(frm)
    rb.grid(row=5, column=0, columnspan=2, sticky="ew")
    runb = ttk.Button(rb, text="Run analysis")
    runb.pack(side="left")
    openb = ttk.Button(rb, text="Open report", state="disabled",
                       command=lambda: state["report"] and webbrowser.open(Path(state["report"]).as_uri()))
    openb.pack(side="left", padx=6)
    ttk.Button(rb, text="Open output folder",
               command=lambda: webbrowser.open(Path(outv.get()).expanduser().as_uri())).pack(side="left")
    pb = ttk.Progressbar(rb, mode="indeterminate", length=200)
    pb.pack(side="right")
    log = tk.Text(frm, height=12, wrap="word", state="disabled")
    log.grid(row=6, column=0, columnspan=2, sticky="nsew", pady=(8, 0))
    frm.rowconfigure(6, weight=1)
    frm.rowconfigure(1, weight=1)

    def write(msg: str):
        log.configure(state="normal")
        log.insert("end", msg + "\n")
        log.see("end")
        log.configure(state="disabled")

    def worker(inputs, outdir, opt, threads, mrg):
        try:
            samples = discover_samples(inputs, merge=mrg)
            q.put(("log", f"{len(samples)} sample(s) found"))
            run(samples, outdir, opt, threads=threads, progress=lambda m: q.put(("log", m)))
            q.put(("done", str(Path(outdir) / "report.html")))
        except Exception as e:  # noqa: BLE001
            q.put(("error", f"{e}\n{traceback.format_exc()}"))

    def start():
        if state["running"]:
            return
        inputs = list(lst.get(0, "end"))
        if not inputs:
            messagebox.showwarning("HuntingTool", "Add at least one FASTQ file or folder.")
            return
        try:
            opt = Options(mode=modev.get(), motif=motifv.get().upper(), min_length=int(minlen.get()),
                          min_q=float(minq.get()), min_fraction=float(minfrac.get()), write_reads=wreads.get())
            threads = int(thr.get())
        except (tk.TclError, ValueError):
            messagebox.showerror("HuntingTool", "Check the parameter values.")
            return
        state["running"] = True
        runb.configure(state="disabled")
        openb.configure(state="disabled")
        pb.start(12)
        write(f"HuntingTool {__version__} – output: {outv.get()}")
        threading.Thread(target=worker, args=(inputs, str(Path(outv.get()).expanduser()), opt, threads, merge.get()),
                         daemon=True).start()

    runb.configure(command=start)

    def poll():
        try:
            while True:
                kind, msg = q.get_nowait()
                if kind == "log":
                    write(msg)
                else:
                    state["running"] = False
                    pb.stop()
                    runb.configure(state="normal")
                    if kind == "done":
                        state["report"] = msg
                        openb.configure(state="normal")
                        write("Done. Click 'Open report' to view the results.")
                    else:
                        write("ERROR: " + msg)
                        messagebox.showerror("HuntingTool", msg.splitlines()[0])
        except queue.Empty:
            pass
        root.after(150, poll)

    for i in inputs or []:
        lst.insert("end", i)
    if outdir:
        outv.set(outdir)
    if autorun:
        root.after(300, start)
    if autoclose_ms:
        root.after(autoclose_ms, root.destroy)
    poll()
    root.mainloop()
    return 0
