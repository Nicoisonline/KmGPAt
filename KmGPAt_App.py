import customtkinter as ctk
from PIL import Image, ImageTk
import ncbi_interactions as ncbi
import threading
import webbrowser
import os
import kmer_analysis as kt

class KmGPAt_App(ctk.CTk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title('KmGPAt')
        self.geometry('1000x600')
        #self.resizable(False, False)

        self.grid_columnconfigure((1), weight=1)
        self.grid_rowconfigure((0), weight=1)

        #Menu
        self.menu = ctk.CTkFrame(self)
        self.menu.grid(row=0, column=0, sticky = "nsew", rowspan=1)
        

        self.menutitle = ctk.CTkLabel(self.menu, text="KmGPAt", font=("Arial", 24))
        self.menutitle.grid(row=0, column=0, sticky = "ew")

        self.bouton_download = ctk.CTkButton(self.menu, text="Download", command=self.open_download)
        self.bouton_download.grid(row=1, column=0, sticky = "ew")

        self.toplevel_download = None

        self.label_kmer = ctk.CTkLabel(self.menu, text="Kmer size: ")
        self.label_kmer.grid(row=5, column=0, sticky = "ew")

        self.show_option = ctk.CTkCheckBox(self.menu, text="Show with Matplotlib", variable=ctk.IntVar(value=0))
        self.show_option.grid(row=4, column=0, sticky = "ew")

        self.entry_kmer = ctk.CTkEntry(self.menu)
        self.entry_kmer.grid(row=6, column=0, sticky = "ew")

        self.bouton_kmer = ctk.CTkButton(self.menu, text="Compute", command=self.compute_kmer)
        self.bouton_kmer.grid(row=7, column=0, sticky = "ew")

        self.files_menu()


        #Window
        self.windows_image = ctk.CTkFrame(self)

        try:
            self.image = ctk.CTkImage(Image.open("output.png"), Image.open("output.png"), size=(578, 417))
            self.label_image = ctk.CTkLabel(self.windows_image, image=self.image, text="")
            self.label_image.grid(row=0, column=0)
        except:
            self.label_image = ctk.CTkLabel(self.windows_image, text="No image available")
            self.label_image.grid(row=0, column=0)

        self.windows_image.grid(row=0, column=1, sticky = "")

        #Data Frame

        self.data_frame()

        self.mainloop()

    def data_frame(self):
        self.data_frame = ctk.CTkFrame(self)
        self.data_frame.grid(row=0, column=2, sticky = "nsew")

        self.data_title = ctk.CTkLabel(self.data_frame, text="Data", font=("Arial", 24))
        self.data_title.grid(row=0, column=0, sticky = "ew")

        self.data_text = ctk.CTkTextbox(self.data_frame, height=300)
        self.data_text.grid(row=1, column=0, sticky = "sew")

    def update_data(self, data):
        self.data_text.delete("0.0", "end")
        self.data_text.insert("0.0", text=data)

    def compute_kmer(self):
        """Compute the kmer and display the image"""
        show_status = bool(self.show_option.get())
        result = "No results available"
        if self.type.get() == "Compare":
            if self.SegBouton_gen_prot.get() == "Genomes":

                result = kt.kmer_pipeline(self.file_1_combobox.get()[:-4], self.file_2_combobox.get()[:-4], int(self.entry_kmer.get()), save=True, comparaison_mode=self.comparaison_mode.get(), show=show_status)
            else:
                result = kt.kmer_pipeline_amino_acids(self.file_1_combobox.get()[:-4], self.file_2_combobox.get()[:-4], int(self.entry_kmer.get()), save=True, comparaison_mode=self.comparaison_mode.get(), show=show_status)
        try:
            self.image = ctk.CTkImage(Image.open("output.png"), Image.open("output.png"), size=(578, 417))
            self.label_image.configure(image=self.image, text="")
        except:
            self.label_image.configure(text="No image available")



        if self.type.get() == "Single":
            if self.window_or_kmer_select.get() == "Window":
                heatmap_or_variance_int = 0
                if self.heatmap_or_variance.get() == "Variance":
                    heatmap_or_variance_int = 1
                # Window
                if self.singleSegBouton_gen_prot_value.get() == "Genomes":
                    result = kt.single_pipeline(self.file_combobox.get()[:-4], int(self.entry_kmer.get()), kmer_or_window=1,save=True, window_size=int(self.window_size_entry.get()), heatmap_or_variance=heatmap_or_variance_int, show=show_status)
                else:
                    result = kt.single_pipeline_amino_acids(self.file_combobox.get()[:-4], int(self.entry_kmer.get()), kmer_or_window=1,save=True, window_size=int(self.window_size_entry.get()), heatmap_or_variance=heatmap_or_variance_int, show=show_status)
            else:
                # Kmer
                if self.singleSegBouton_gen_prot_value.get() == "Genomes":
                    result = kt.single_pipeline(self.file_combobox.get()[:-4], int(self.entry_kmer.get()), kmer_or_window=0,save=True, window_size=None, heatmap_or_variance=None, show=show_status)
                else:
                    result = kt.single_pipeline_amino_acids(self.file_combobox.get()[:-4], int(self.entry_kmer.get()), kmer_or_window=0,save=True, window_size=None, heatmap_or_variance=None, show=show_status)
            try:
                self.image = ctk.CTkImage(Image.open("output.png"), Image.open("output.png"), size=(578, 417))
                self.label_image.configure(image=self.image, text="")
            except:
                self.label_image.configure(text="No image available")
        self.update_data(result)

    def files_menu(self):

        # Variables for names
        self.genomes_list = os.listdir("data/genomes")
        self.protseq_list = os.listdir("data/protseq")
        self.genomes_list_show = self.genomes_list.copy()
        self.protseq_list_show = self.protseq_list.copy()

        # Type

        self.type = ctk.StringVar(value="Compare")
        self.type_combobox = ctk.CTkSegmentedButton(self.menu, values=["Compare", "Single"], variable=self.type, command=self.frame_switch)
        self.type_combobox.grid(row=2, column=0, sticky = "ew")

        # Compare

        self.compare_frame = ctk.CTkFrame(self.menu)
        self.compare_frame.grid(row=3, column=0, sticky = "ew")

        self.SegBouton_gen_prot_value = ctk.StringVar(value="Genomes")
        self.SegBouton_gen_prot = ctk.CTkSegmentedButton(self.compare_frame, values=["Genomes", "ProtSeq"], variable=self.SegBouton_gen_prot_value, command=self.combobox_set)
        self.SegBouton_gen_prot.grid(row=6, column=0, sticky = "")

        self.searchbar_title = ctk.CTkLabel(self.compare_frame, text="Search: ")
        self.searchbar_title.grid(row=7, column=0, sticky = "ew")

        self.researchbar = ctk.CTkEntry(self.compare_frame)
        self.researchbar.bind("<KeyRelease>", self.combobox_set)
        self.researchbar.grid(row=8, column=0, sticky = "ew")

        self.file_1_title = ctk.CTkLabel(self.compare_frame, text="File 1: ")
        self.file_1_title.grid(row=9, column=0, sticky = "ew")

        self.file_1_combobox = ctk.CTkOptionMenu(self.compare_frame, values=self.genomes_list)
        self.file_1_combobox.grid(row=10, column=0, sticky = "")
        #self.file_1_combobox.bind("<Enter>", lambda event: self.file_1_combobox.configure(values=os.listdir("data/genomes")))

        self.file_2_title = ctk.CTkLabel(self.compare_frame, text="File 2: ")
        self.file_2_title.grid(row=11, column=0, sticky = "ew")

        self.file_2_combobox = ctk.CTkOptionMenu(self.compare_frame, values=self.genomes_list)
        self.file_2_combobox.grid(row=12, column=0, sticky = "")
        #self.file_2_combobox.bind("<Enter>", lambda event: self.file_2_combobox.configure(values=os.listdir("data/genomes")))

        self.comparaison_mode = ctk.CTkCheckBox(self.compare_frame, text="Switch Mode", variable=ctk.IntVar(value=0))
        self.comparaison_mode.grid(row=13, column=0, sticky = "ew")

        # Single

        self.single_frame = ctk.CTkFrame(self.menu)
        self.single_frame.grid(row=3, column=0, sticky = "ew")

        self.singleSegBouton_gen_prot_value = ctk.StringVar(value="Genomes")
        self.singleSegBouton_gen_prot = ctk.CTkSegmentedButton(self.single_frame, values=["Genomes", "ProtSeq"], variable=self.singleSegBouton_gen_prot_value, command=self.single_combobox_set)
        self.singleSegBouton_gen_prot.grid(row=6, column=0, sticky = "")

        self.single_searchbar_title = ctk.CTkLabel(self.single_frame, text="Search: ")
        self.single_searchbar_title.grid(row=7, column=0, sticky = "ew")

        self.single_researchbar = ctk.CTkEntry(self.single_frame)
        self.single_researchbar.bind("<KeyRelease>", self.single_combobox_set)
        self.single_researchbar.grid(row=8, column=0, sticky = "ew")

        self.file_title = ctk.CTkLabel(self.single_frame, text="File: ")
        self.file_title.grid(row=9, column=0, sticky = "ew")

        self.file_combobox = ctk.CTkOptionMenu(self.single_frame, values=self.genomes_list)
        self.file_combobox.grid(row=10, column=0, sticky = "")

        self.window_or_kmer_select = ctk.CTkSegmentedButton(self.single_frame, values=["Kmer", "Window"], variable=ctk.StringVar(value="Kmer"), command=self.window_or_kmer)
        self.window_or_kmer_select.grid(row=11, column=0, sticky = "ew")

        self.window_size = ctk.CTkLabel(self.single_frame, text="Window size: ")

        self.window_size_entry = ctk.CTkEntry(self.single_frame)

        self.heatmap_or_variance = ctk.CTkSegmentedButton(self.single_frame, values=["Heatmap", "Variance"], variable=ctk.StringVar(value="Heatmap"))

        self.save_window_label = ctk.CTkLabel(self.single_frame, text="Save window", font=("Arial", 12, "bold"))

        self.window_number_label = ctk.CTkLabel(self.single_frame, text="Window number: ")

        self.window_number_entry = ctk.CTkEntry(self.single_frame)

        self.save_button = ctk.CTkButton(self.single_frame, text="Save", command=self.save_window)

        self.PCA_button = ctk.CTkButton(self.single_frame, text="RunPCA", command=self.compute_PCA)

        self.PCA_disclaimer = ctk.CTkLabel(self.single_frame, text="PCA - No Protseq", font=("Arial", 12, "bold"))

        self.single_frame.grid_forget()

    def compute_PCA(self):
        kt.PCA_pipeline(self.file_combobox.get()[:-4], k=int(self.entry_kmer.get()), window_number=int(self.window_size_entry.get()), genomes_or_protseq=self.singleSegBouton_gen_prot.get(), show=bool(self.show_option.get()), save=True)
        # We update the image
        try:
            self.image = ctk.CTkImage(Image.open("output.png"), Image.open("output.png"), size=(578, 417))
            self.label_image.configure(image=self.image, text="")
        except:
            self.label_image.configure(text="No image available")

    def save_window(self):
        # Si génome :
        if self.singleSegBouton_gen_prot_value.get() == "Genomes":
            kt.save_windows(self.file_combobox.get()[:-4], int(self.window_size_entry.get()), int(self.window_number_entry.get()), genome_or_protseq=0)
        else:
            kt.save_windows(self.file_combobox.get()[:-4], int(self.window_size_entry.get()), int(self.window_number_entry.get()), genome_or_protseq=1)

    def window_or_kmer(self, event):
        if self.window_or_kmer_select.get() == "Window":
            self.window_size.grid(row=12, column=0, sticky = "ew")
            self.window_size_entry.grid(row=13, column=0, sticky = "ew")
            self.heatmap_or_variance.grid(row=14, column=0, sticky = "ew")
            self.save_window_label.grid(row=15, column=0, sticky = "ew")
            self.window_number_label.grid(row=16, column=0, sticky = "ew")
            self.window_number_entry.grid(row=17, column=0, sticky = "ew")
            self.save_button.grid(row=18, column=0, sticky = "ew")
            self.PCA_button.grid(row=19, column=0, sticky = "ew")
            self.PCA_disclaimer.grid(row=20, column=0, sticky = "ew")
        else:
            # Kmer
            self.window_size.grid_forget()
            self.window_size_entry.grid_forget()
            self.heatmap_or_variance.grid_forget()
            self.save_window_label.grid_forget()
            self.window_number_label.grid_forget()
            self.window_number_entry.grid_forget()
            self.save_button.grid_forget()
            self.PCA_button.grid_forget()
            self.PCA_disclaimer.grid_forget()


    def update_list(self):
        self.genomes_list = os.listdir("data/genomes")
        self.protseq_list = os.listdir("data/protseq")
        self.genomes_list_show = self.genomes_list.copy()
        self.protseq_list_show = self.protseq_list.copy()

    def frame_switch(self, event):
        if self.type.get() == "Compare":
            self.compare_frame.grid(row=3, column=0, sticky = "sw")
            self.single_frame.grid_forget()
            self.combobox_set(None)
        else:
            self.single_frame.grid(row=3, column=0, sticky = "sw")
            self.compare_frame.grid_forget()
            self.single_combobox_set(None)

    def single_combobox_set(self, event):
        self.update_list()
        # We check if the research bar is empty
        if self.single_researchbar.get() == "":
            if self.singleSegBouton_gen_prot_value.get() == "Genomes":
                self.file_combobox.configure(values=self.genomes_list)
                
            else:
                self.file_combobox.configure(values=self.protseq_list)

        # We check if the research bar is not empty
        else:
            if self.singleSegBouton_gen_prot_value.get() == "Genomes":
                self.genomes_list_show = [x for x in self.genomes_list if self.single_researchbar.get() in x]
                self.file_combobox.configure(values=self.genomes_list_show)

            else:
                self.protseq_list_show = [x for x in self.protseq_list if self.single_researchbar.get() in x]
                self.file_combobox.configure(values=self.protseq_list_show)

    def combobox_set(self, event):
        self.update_list()
        # We check if the research bar is empty
        if self.researchbar.get() == "":
            if self.SegBouton_gen_prot_value.get() == "Genomes":
                self.file_1_combobox.configure(values=self.genomes_list)
                self.file_2_combobox.configure(values=self.genomes_list)
                
            else:
                self.file_1_combobox.configure(values=self.protseq_list)
                self.file_2_combobox.configure(values=self.protseq_list)

        # We check if the research bar is not empty
        else:
            if self.SegBouton_gen_prot_value.get() == "Genomes":
                self.genomes_list_show = [x for x in self.genomes_list if self.researchbar.get() in x]
                self.file_1_combobox.configure(values=self.genomes_list_show)
                self.file_2_combobox.configure(values=self.genomes_list_show)

            else:
                self.protseq_list_show = [x for x in self.protseq_list if self.researchbar.get() in x]
                self.file_1_combobox.configure(values=self.protseq_list_show)
                self.file_2_combobox.configure(values=self.protseq_list_show)

    def open_download(self):
        self.update_list()
        if self.toplevel_download is None or not self.toplevel_download.winfo_exists():
            self.toplevel_download = KmGPAt_Download(self)
            self.toplevel_download.after(200, self.toplevel_download.lift)
        else:
            self.toplevel_download.focus()



class KmGPAt_Download(ctk.CTkToplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title('KmGPAt - Download')
        self.geometry('400x300')
        self.resizable(False, False)

        self.grid_columnconfigure((0), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)

        self.title = ctk.CTkLabel(self, text="Download", font=("Arial", 24))
        self.title.grid(row=0, column=0, sticky = "nesw")

        self.text = ctk.CTkLabel(self, text="Please provide the TaxName or Accession")
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
        self.progressbar.start()
        thread = threading.Thread(target=self.thread_dl)
        thread.start()

    def thread_dl(self):
        if self.input_dl.get() == "summary":
            res = ncbi.download_summary_bacteria()
        else:
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
    ncbi.check_files_structure()
    ctk.set_appearance_mode("dark")
    KmGPAt_App()

if __name__ == "__main__":
    main()