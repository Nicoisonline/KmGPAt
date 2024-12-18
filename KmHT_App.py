import customtkinter as ctk
import ncbi_interactions as ncbi
import threading
import webbrowser

class KmHT_App(ctk.CTk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title('KmHT')
        self.geometry('800x600')
        self.resizable(False, False)

        self.grid_columnconfigure((1), weight=1)

        #Menu
        self.menutitle = ctk.CTkLabel(self, text="KmHT", font=("Arial", 24))
        self.menutitle.grid(row=0, column=0, sticky = "sw")

        self.bouton_download = ctk.CTkButton(self, text="Download", command=self.open_download)
        self.bouton_download.grid(row=1, column=0, sticky = "sw")

        self.toplevel_download = None

        self.mainloop()

    def open_download(self):
        if self.toplevel_download is None or not self.toplevel_download.winfo_exists():
            self.toplevel_download = KmHT_Download(self)
            self.toplevel_download.after(200, self.toplevel_download.lift)
        else:
            self.toplevel_download.focus()



class KmHT_Download(ctk.CTkToplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title('KmHT - Download')
        self.geometry('400x300')
        self.resizable(False, False)

        self.grid_columnconfigure((0), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)

        self.title = ctk.CTkLabel(self, text="Download", font=("Arial", 24))
        self.title.grid(row=0, column=0, sticky = "nesw")

        self.text = ctk.CTkLabel(self, text="Please provide the taxname or projectID")
        self.text.grid(row=1, column=0, sticky = "ew")

        self.entry = ctk.CTkEntry(self)
        self.entry.grid(row=2, column=0, sticky = "ew")

        self.bouton = ctk.CTkButton(self, text="Download", command=self.download)
        self.bouton.grid(row=3, column=0, sticky = "")

        self.status = ctk.CTkLabel(self, text="Status: Ready", bg_color="grey", corner_radius=5)
        self.status.grid(row=4, column=0, sticky = "ew")

        self.progressbar = ctk.CTkProgressBar(self, orientation="horizontal", mode="indeterminate")
        self.progressbar.grid(row=5, column=0, sticky = "ew")

        self.ncbi_link = "https://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/"

        self.source = ctk.CTkLabel(self, text="Source: National Center for Biotechnology Information (NCBI).\nArchive of Old RefSeq Bacterial Genomes.")
        self.source.grid(row=6, column=0, sticky = "ew")

        self.source._label.bind("<Button-1>", lambda event: webbrowser.open_new_tab(self.ncbi_link))
        self.source._label.bind("<Enter>", lambda event: self.source.configure(font=ctk.CTkFont(underline=True), cursor="hand2"))
        self.source._label.bind("<Leave>", lambda event: self.source.configure(font=ctk.CTkFont(underline=False), cursor="arrow"))

    def download(self):
        if hasattr(self, 'downloading') and self.downloading:
            self.status.configure(text="Status: Download already in progress")
            return

        self.downloading = True
        self.input_dl = ctk.StringVar(value=self.entry.get())
        self.status.configure(text="Status: Downloading - This can take a while")
        if self.input_dl.get() == "summary":
            print("Downloading bacteria summary")
            ncbi.download_summary_bacteria()
            self.downloading = False
        else:
            print("Downloading bacteria data")
            self.progressbar.start()
            thread = threading.Thread(target=self.thread_dl)
            thread.start()

    def thread_dl(self):
        res = ncbi.download_bacteria(str(self.input_dl.get()))
        if res == -1:
            self.status.configure(text="Status: Please provide the taxname or projectID")
        if res == 0:
            self.status.configure(text="Status: Download completed")
        if res == 1:
            self.status.configure(text="Status: Protfile already downloaded, genomes file downloaded")
        if res == 2:
            self.status.configure(text="Status: Protfile downloaded, genomes file downloaded")
        if res == 3 or res == 4:
            self.status.configure(text="Status: Error downloading file")
        if res == 12:
            self.status.configure(text="Status: All files already downloaded")
        self.progressbar.stop()
        self.downloading = False

def main():
    ctk.set_appearance_mode("dark")
    KmHT_App()

if __name__ == "__main__":
    main()