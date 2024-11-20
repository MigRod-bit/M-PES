import tkinter as tk
import customtkinter
import csv
from tkinter import messagebox, ttk, filedialog

def show_selection():
    try:
        # Obtener el ID del primer elemento seleccionado.
        item = treeview.selection()[0]
    except IndexError:
        # Si la tupla está vacía, entonces no hay ningún
        # elemento seleccionado.
        messagebox.showwarning(
            message="Debe seleccionar un elemento.",
            title="No hay selección"
        )
    else:
        # A partir de este ID retornar el texto del elemento.
        text = treeview.item(item, option="text")
        # Mostrarlo en un cuadro de diálogo.
        messagebox.showinfo(message=text, title="Selección")

main_window = tk.Tk()
main_window.title("Vista de árbol en Tkinter")
treeview = ttk.Treeview(columns=("size", "lastmod"))
treeview.heading("#0", text="Archivo")
treeview.heading("size", text="Tamaño")
treeview.heading("lastmod", text="Última modificación")
treeview.insert(
    "", 
    tk.END,
    text="Hola.exe", 
    values=("76 Bt", "7.30")
)
treeview.insert(
    "", 
    tk.END,
    text="adios.exe", 
    values=("796 Bt", "17.30")
)
treeview.pack()
button = ttk.Button(text="Mostrar selección",
                    command=show_selection)
button.pack()
main_window.mainloop()

