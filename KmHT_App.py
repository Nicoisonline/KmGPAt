import customtkinter as ctk
import ncbi_interactions as ncbi
import time
import threading

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
        self.title.grid(row=0, column=0, sticky = "sw")

        self.text = ctk.CTkLabel(self, text="Please provide the taxname or projectID")
        self.text.grid(row=1, column=0, sticky = "sw")

        self.entry = ctk.CTkEntry(self)
        self.entry.grid(row=2, column=0, sticky = "sw")

        self.bouton = ctk.CTkButton(self, text="Download", command=self.download)
        self.bouton.grid(row=3, column=0, sticky = "sw")

        self.status = ctk.CTkLabel(self, text="Status: Ready")
        self.status.grid(row=4, column=0, sticky = "sw")

        self.progressbar = ctk.CTkProgressBar(self, orientation="horizontal", mode="indeterminate")
        self.progressbar.grid(row=5, column=0, sticky = "sw")

    def download(self):
        input = self.entry.get()
        self.status.configure(text="Status: Downloading")
        self.progressbar.start()
        if input == "summary":
            print("Downloading bacteria summary")
            ncbi.download_summary_bacteria()
        else:
            print("Downloading bacteria data")
            threading.Thread(target=self.download_func(input)).start()

        
    def download_func(self, input):
        ncbi.download_bacteria(input)
        self.progressbar.stop()
        self.status.configure(text="Status: Downloaded Success - Ready")

        


def main():
    ctk.set_appearance_mode("dark")
    KmHT_App()

if __name__ == "__main__":
    main()