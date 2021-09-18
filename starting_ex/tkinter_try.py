from tkinter import *

#create Root widget which is the window
root = Tk()
#Simple label example with .pack()
'''
#Create label widget 
my_label = Label(root, text='Hello World!')
#my_label.delete() #to delete text
#Append label to the center of the window 
my_label.pack()
'''

#Simple label example with .grid()
'''
#Create widget
my_label1 = Label(root, text='Hello World!')
my_label2 = Label(root, text='How u doing').grid(row=1,column=1)
#Append on the screen
my_label1.grid(row=0, column=0)#columnspan=int (in case multiple columns)
'''

#Create button widget and implement function LAMBDA EXAMPLE 
'''
def my_click():
    my_label = Label(root, text='You clicked!').pack()
my_button = Button(root, text='Click me!', command=my_click) #state=DISABLED, padx/pady=int for size,
#in case you want to send parameters
lambda: my_click(0,2)
my_button.pack()
'''

#Create input box (Entry object)
'''
e = Entry(root) #width=int, bg='blue'(background), fg='colour'(foreground),
e.pack()
#Give default text
e.insert(0, 'Enter your name')
#e.get() gets whatever is written in the box
#Used in a function attached to a button
def my_click():
    text = 'Hello ' + e.get()
    my_label = Label(root, text=text).pack()
#Append functional button
my_button = Button(root, text='Say hello!', command=my_click).pack()
'''


#Add title and icon to program
'''
root.title('My program title')
#for icon check accepted files
#root.iconbitmap('/path/to/file')
'''

#Using images! Need to import a module to manage different formats
'''
#the module is pillow but it's referenced ad PIL
from PIL import Image, ImageTk
#Open the image, pass it to photoImage object 
my_image = ImageTk.PhotoImage(Image.open('computer_icon.png'))
#Add image to label and append it
image_label = Label(image=my_image).pack()
'''







#Create event loop to start program
#root.mainloop()

query_length = 300
difference_from_query = 50
lower, upper = query_length*(difference_from_query/100), query_length*((difference_from_query + 100)/100)
print(lower, upper)