import pyglet
pyglet.options['debug_gl'] = False
import math
import numpy as np
import tkinter as tk
from tkinter import filedialog as fd
import loadingroutine as lr
from collections import deque
import multiprocessing as mul
from tkinter import messagebox as errnotify
import sys
sys.path.append('../')
from herdingspikes.probe import NeuroPixel, BioCam
import copy
import time
import random
import gc
import os


if __name__ == "__main__":

    post_processed = False
    SETTINGS_MODE = False

    if not os.path.isdir("gen"):
        os.mkdir('gen')

    def center(win, win2):
        win.update_idletasks()
        win2.update_idletasks()
        print(win2.winfo_width(), win.winfo_width())
        width = win2.winfo_reqwidth()
        frm_width = win.winfo_rootx() - win.winfo_x()
        win_width = width + 2 * frm_width
        height = win.winfo_height()
        titlebar_height = win.winfo_rooty() - win.winfo_y()
        win_height = height + titlebar_height + frm_width
        x = win.winfo_screenwidth() // 2 - win_width // 2
        y = win.winfo_screenheight() // 2 - win_height // 2
        win.geometry('{}x{}+{}+{}'.format(width, height, x, y))
        win.deiconify()

    STEPS_TO_SEE_MAX = 400
    STEPS_TO_SEE_MAX = max(1, STEPS_TO_SEE_MAX)

    THEME = 0  # 0: Dark, 1: Bright

    file_raw = None
    file_proc = None

    folder_name = None

    is_biocam = False

    skip_loading = False

    root = tk.Tk()
    img = tk.PhotoImage(file='icon.gif')
    root.tk.call('wm', 'iconphoto', root._w, img)
    root.protocol("WM_DELETE_WINDOW", lambda: exit())

    class TkSetup:

        def __init__(self, parent):

            self.parent = parent
            self.labels = []
            self.widgets = []

            self.probe_mode_var = tk.StringVar(self.parent, value="Neuropixels")
            self.labels.append(tk.Label(self.parent, text='Probe:'))
            self.widgets.append(tk.OptionMenu(self.parent, self.probe_mode_var, 'Neuropixels', 'Biocam'))

            self.labels.append(tk.Label(self.parent, text='STS_MAX:'))
            self.widgets.append(tk.Frame(self.parent))
            n_frame = self.widgets[len(self.widgets) - 1]
            self.STS_MAX_SET = 100
            self.STS_MAX_text = tk.Text(n_frame, width=16, height=1)
            self.STS_MAX_text.grid(row=0, column=0)
            self.STS_MAX_text.insert(tk.INSERT, str(self.STS_MAX_SET))

            self.labels.append(tk.Label(self.parent, text='* Raw Signal File:'))
            self.widgets.append(tk.Frame(self.parent))
            n_frame_2 = self.widgets[len(self.widgets) - 1]
            self.raw_file = None
            self.raw_file_button = tk.Button(n_frame_2, text="Open...", command=self.open_raw)
            self.raw_file_button.grid(row=0, column=0)
            self.raw_file_text = tk.Label(n_frame_2, text="")
            self.raw_file_text.grid(row=0, column=1)

            self.labels.append(tk.Label(self.parent, text='* Processed Spikes File:'))
            self.widgets.append(tk.Frame(self.parent))
            n_frame_3 = self.widgets[len(self.widgets) - 1]
            self.proc_file = None
            self.proc_file_button = tk.Button(n_frame_3, text="Open...", command=self.open_proc)
            self.proc_file_button.grid(row=0, column=0)
            self.proc_file_text = tk.Label(n_frame_3, text="")
            self.proc_file_text.grid(row=0, column=1)

            self.labels.append(tk.Label(self.parent, text='Preprocessed SpikeSeer Folder:'))
            self.widgets.append(tk.Frame(self.parent))
            n_frame_4 = self.widgets[len(self.widgets) - 1]
            self.proc_folder = None
            self.proc_folder_button = tk.Button(n_frame_4, text="Open...", command=self.open_fldr)
            self.proc_folder_button.grid(row=0, column=0)
            self.proc_folder_text = tk.Label(n_frame_4, text="")
            self.proc_folder_text.grid(row=0, column=1)

            self.gfx_mode_var = tk.StringVar(self.parent, value="Pyglet")
            self.labels.append(tk.Label(self.parent, text='Graphics Framework:'))
            self.widgets.append(tk.OptionMenu(self.parent, self.gfx_mode_var, 'Pyglet'))

            self.theme_var = tk.StringVar(self.parent, value="Dark")
            self.labels.append(tk.Label(self.parent, text='Theme (Pyglet):'))
            self.widgets.append(tk.OptionMenu(self.parent, self.theme_var, 'Dark'))

            for cnt, label in enumerate(self.labels):
                label.grid(column=0, row=cnt, sticky=tk.E, padx=4)

            for cnt, widget in enumerate(self.widgets):
                widget.grid(column=1, row=cnt)

            self.OK = tk.Button(parent, text="Confirm", command=self.confirm)
            self.OK.grid(column=0, row=cnt+1, columnspan=2, pady=10)

        def open_raw(self):
            self.raw_file = fd.askopenfilename(title="Open Raw Signal File", initialdir='./')
            if self.raw_file is not None:
                self.raw_file_text.config(text=self.raw_file.replace('\\', '/').split('/')[-1])
            center(root, self.parent)

        def open_proc(self):
            self.proc_file = fd.askopenfilename(title="Open Processed Spikes File", initialdir='./')
            if self.proc_file is not None:
                self.proc_file_text.config(text=self.proc_file.replace('\\', '/').split('/')[-1])
            center(root, self.parent)

        def open_fldr(self):
            self.proc_folder = fd.askdirectory(title="Open Generated SpikeSeer Folder", initialdir='./gen/')
            if self.proc_folder is not None:
                self.proc_folder_text.config(text=self.proc_folder.replace('\\', '/').split('/')[-1])
            center(root, self.parent)

        def confirm(self):
            global skip_loading, STEPS_TO_SEE_MAX, file_raw, file_proc, THEME, is_biocam, folder_name
            prb = self.probe_mode_var.get()
            if prb == 'Biocam':
                is_biocam = True
            sts = 1
            try:
                sts = int(self.STS_MAX_text.get("1.0", tk.END).strip())
            except Exception:
                errnotify.showerror("SpikeSeer: Conversion Error",
                                    "STS_MAX could not be interpreted as an integer. Aborting.")
                exit()
            STEPS_TO_SEE_MAX = sts
            file_raw = self.raw_file
            file_proc = self.proc_file
            proc_fldr = self.proc_folder
            gfx_framework = self.gfx_mode_var.get()
            theme = self.theme_var.get()
            if theme == 'Bright':
                THEME = 1
            if proc_fldr is None or proc_fldr == '':
                self.t_lvl = tk.Toplevel()
                self.t_lvl.protocol("WM_DELETE_WINDOW", lambda: exit())
                img = tk.PhotoImage(file='icon.gif')
                self.t_lvl.tk.call('wm', 'iconphoto', self.t_lvl._w, img)
                _ = tk.Label(self.t_lvl, text='Input folder name for generated files')
                _.pack()
                self.txt = tk.Text(self.t_lvl, height=1, width=30)
                self.txt.pack()
                btn = tk.Button(self.t_lvl, text='Confirm', command=self.apl_proc_fldr)
                btn.pack()
                center(self.t_lvl, self.t_lvl)
            else:
                skip_loading = True
                folder_name = proc_fldr
                root.destroy()

        def apl_proc_fldr(self):
            global folder_name
            folder_name = self.txt.get("1.0", tk.END)
            self.t_lvl.destroy()
            root.destroy()

    frm = tk.Frame(root)
    root.wm_title("Setup")
    _ = TkSetup(frm)
    frm.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    center(root, _.parent)
    root.mainloop()

    def gen_smoothed_arrays(a):

        _scalar = 3000

        a0 = (a / np.amax(a)) * _scalar

        a1 = np.log(a + 1)
        a1 = (a1 / np.amax(a1)) * _scalar

        a2 = np.log(np.log(a + 1) + 1)
        a2 = (a2 / np.amax(a2)) * _scalar

        a3 = np.log(np.log(np.log(a + 1) + 1) + 1)
        a3 = (a3 / np.amax(a3)) * _scalar

        a4 = np.log(np.log(np.log(np.log(a + 1) + 1) + 1) + 1)
        a4 = (a4 / np.amax(a4)) * _scalar

        a5 = np.log(np.log(np.log(np.log(np.log(a + 1) + 1) + 1) + 1) + 1)
        a5 = (a5 / np.amax(a5)) * _scalar

        return [a0, a1, a2, a3, a4, a5]

    brightness_bank = []
    for color in range(0, 1000 * 3 + 1):
        r = 0
        g = 0
        b = 0
        bank = color/1000
        r = min(1., bank)
        bank -= r
        g = min(1., bank)
        bank -= g
        b = bank
        brightness_bank.append([r, g, b])

    brightness_bank_2 = []
    for color in range(0, 1000 * 3 + 1):
        r = 0
        g = 0
        b = 0
        a = 0
        if 0 <= color < 1000:
            b = 0.5
            a = color/1000
        if 1000 <= color < 2000:
            r = (color - 1000)/1000
            b = 0.5
            a = 1
        if 2000 <= color <= 3000:
            g = (color - 2000)/1000
            r = 1
            b = 0.5
            a = 1

        brightness_bank_2.append([r, g, b, a])

    b_THRESH = 1
    b_DECAY = 24

    class SpikeSquare:

        def __init__(self, x_, y_, w, h):
            self.x = x_
            self.y = y_
            self.w = w
            self.h = h
            self.brightness = 0
            self.m_brightness = 1000 * 3

        def draw(self):

            if self.brightness < b_THRESH:
                return

            pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
            bank_item = brightness_bank[int(self.brightness)]
            r, g, b = bank_item[0], bank_item[1], bank_item[2]
            pyglet.gl.glColor3f(r, g, b)

            pyglet.gl.glVertex2f(half_w + int((self.x + ox) * zoom), half_h + int((self.y + oy) * zoom))
            pyglet.gl.glVertex2f(half_w + int((self.x + ox + self.w) * zoom), half_h + int((self.y + oy) * zoom))
            pyglet.gl.glVertex2f(half_w + int((self.x + ox + self.w) * zoom), half_h + int((self.y + oy + self.h) * zoom))
            pyglet.gl.glVertex2f(half_w + int((self.x + ox) * zoom), half_h + int((self.y + oy + self.h) * zoom))

            pyglet.gl.glEnd()

        def update(self):
            self.brightness -= b_DECAY
            self.brightness = max(0, self.brightness)

        def trigger(self):
            self.brightness = self.m_brightness

        def load(self):
            pass

        def unload(self):
            pass

    def patch_idle_loop():
        def idle(self):
            self.clock.call_scheduled_functions(self.clock.update_time())
            return self.clock.get_sleep_time(True)

        if pyglet.app.EventLoop.idle != idle:
            pyglet.app.EventLoop.idle = idle

    # Window patch for OpenGL
    def patch_window_for_opengl_core():
        def draw_mouse_cursor(self):
            pass

        pyglet.window.BaseWindow.draw_mouse_cursor = draw_mouse_cursor

    patch_idle_loop()
    patch_window_for_opengl_core()

    GLOBAL_SCALAR = 10
    STEPS_TO_SEE = STEPS_TO_SEE_MAX//2
    line_width = 2
    event_index = 0

    root = tk.Tk()
    #root.withdraw()
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    root.geometry("0x0+" + str(int(screen_width*2)) + "+" + str(int(screen_height*2)))

    probe = None

    B_WIDTH = 100
    B_HEIGHT = 100

    if folder_name is None or folder_name.strip() == '':
        errnotify.showerror("SpikeSeer: Missing Folder Name Error",
                            "You either closed the folder name dialog or did not input anything. Aborting program.")
        exit()

    if file_proc is None or file_raw is None or file_proc == '' or file_raw == '':
        errnotify.showerror("SpikeSeer: Missing File Error", "Either raw signal file or processed spikes file missing. Aborting program.")
        exit()

    try:
        print('isbiocam:',is_biocam)
        if is_biocam:
            probe = BioCam(data_file_path=file_raw, fps=8000, masked_channels=[])
            GLOBAL_SCALAR = 200
        else:
            print('is npx')
            probe = NeuroPixel(data_file_path=file_raw, fps=30000, masked_channels=[])
    except AssertionError:
        errnotify.showerror('Probe Type Error', 'Selected probe type threw an exception during instantiation with given raw data file (' + file_raw + '), attempting to launch SpikeSeer with different probe type.')
        try:
            if not is_biocam:
                is_biocam = True
                probe = BioCam(data_file_path=file_raw, fps=8000, masked_channels=[])
                GLOBAL_SCALAR = 200
            else:
                is_biocam = False
                probe = NeuroPixel(data_file_path=file_raw, fps=30000, masked_channels=[])
            errnotify.showinfo('Probe Type Error', 'Different probe model was able to be instantiated with provided raw data file. Proceeding with application.')
        except AssertionError:
            errnotify.showerror('Probe Type Error', 'No probe types successfully instantiated with given raw data. Aborting program.')
            exit()
    data_base = probe.Read(0, 1)
    print(data_base.shape)
    N_FRAMES = probe.nFrames
    MOD_POS = copy.deepcopy(probe.positions)

    channels_linear = []

    min_w = 0
    max_w = 0
    min_h = 0
    max_h = 0

    for elem in MOD_POS:
        elem[0] = elem[0] * GLOBAL_SCALAR
        elem[1] = elem[1] * GLOBAL_SCALAR

    for pos in MOD_POS:
        if pos[0] < min_w:
            min_w = pos[0]
        if pos[0] > max_w:
            max_w = pos[0]

        if pos[1] < min_h:
            min_h = pos[1]
        if pos[1] > max_h:
            max_h = pos[1]

    min_w = min_w // B_HEIGHT
    max_w = max_w // B_HEIGHT
    min_h = min_h // B_WIDTH
    max_h = max_h // B_WIDTH

    # -------------------------------------------------------------------------------------------------------------------

    # Width and height of cell array field
    WIDTH = max_w
    HEIGHT = max_h

    # Contains all the info about the actual signal
    channels = []
    c_gfx = []
    squares = []

    for jnv in range(0, WIDTH + 2):
        channels.append([])
        c_gfx.append([])
        squares.append([])
        for jnh in range(0, HEIGHT + 2):
            squares[jnv].append(SpikeSquare(jnv * B_WIDTH, jnh * B_HEIGHT, B_WIDTH, B_HEIGHT))
            channels[jnv].append([])
            c_gfx[jnv].append([])

    squares_processed = np.zeros((len(squares), len(squares[0])))

    load_end_flag = mul.Value('f',  -1)

    shape1 = mul.Value('i', 0)
    shape2 = mul.Value('i', 0)
    shape3 = mul.Value('i', 0)
    avg = mul.Value('i', 0)

    #COLORS
    BASE_R = 1
    BASE_G = 1
    BASE_B = 1

    BG_R = 0
    BG_G = 0
    BG_B = 0

    if THEME == 1:
        BASE_R = 0
        BASE_G = 0
        BASE_B = 0

        BG_R = 1
        BG_G = 1
        BG_B = 1

        brightness_bank = list(reversed(brightness_bank))

        brightness_bank_2 = []
        for color in range(0, 1000 * 3 + 1):
            r = 0
            g = 0
            b = 0
            a = 0
            if 0 <= color < 1000:
                r = 1
                a = (color/1000)
                b = 1
                g = 0.5
            if 1000 <= color < 1500:
                g = 0.5 - ((color - 1000)/500)
                r = 1
                a = 1
                b = 1
            if 1500 <= color < 2000:
                g = 0
                r = 1 - ((color - 1500)/500)
                a = 1
                b = 1
            if 2000 <= color <= 2500:
                g = 0
                r = 0.5 - ((color - 2000)/500)
                a = 1
                b = 1
            if 2500 <= color <= 3000:
                g = 0
                r = 0
                a = 1
                b = 1 - ((color - 2500)/500)

            brightness_bank_2.append([r, g, b, a])

    #endCOLORS

    pyglet.gl.glClearColor(BG_R, BG_G, BG_B, 1)

    window_radius = 1

    def channel_wipe():
        for chn in channels_linear:
            chn.to_light.clear()
            chn.events.clear()

    def channel_relight():
        for chn in channels_linear:
            chn.relight()

    def log_search():
        global event_index
        drift = 10
        n_thresh = 100
        print("LOG_SEARCH", shape3.value, event_index)
        event_index = min(shape3.value - 1, event_index)
        SPIKE_STEP_TARGET = event_index
        if process_post_event(event_index, mmap_indices, mmap_xy, mmap_data).min_frame_v < CURRENT_FRAME + drift:
            # if abs(process_post_event(self.event_index).min_frame_v - self.STEP - 10) < 100:
            #     return
            jump_step = (shape3.value - event_index) // 2
            sign = 1
            run_flag = True
            while run_flag:
                SPIKE_STEP_TARGET += jump_step * sign
                print(SPIKE_STEP_TARGET)
                new_diff = CURRENT_FRAME + drift - process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v
                new_diff_abs = abs(process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v - CURRENT_FRAME - drift)
                if new_diff_abs < n_thresh or jump_step <= 1:
                    if new_diff < 0:
                        while (CURRENT_FRAME + drift - process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v) < 0:
                            SPIKE_STEP_TARGET -= 1
                        SPIKE_STEP_TARGET += 1
                    else:
                        while (CURRENT_FRAME + drift - process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v) >= 0:
                            SPIKE_STEP_TARGET += 1
                    run_flag = False
                else:
                    print(jump_step, sign)
                    if jump_step <= 1:
                        run_flag = False
                    jump_step = jump_step//2
                    if new_diff < 0:
                        sign = -1
                    else:
                        sign = 1
            event_index = SPIKE_STEP_TARGET

        elif process_post_event(event_index, mmap_indices, mmap_xy, mmap_data).min_frame_v == CURRENT_FRAME + drift:
            pass

        else:
            jump_step = (event_index // 2)
            sign = -1
            run_flag = True
            while run_flag:
                SPIKE_STEP_TARGET += jump_step * sign
                print(SPIKE_STEP_TARGET)

                new_diff = CURRENT_FRAME + drift - process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v
                new_diff_abs = abs(process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v - CURRENT_FRAME - drift)
                if new_diff_abs < n_thresh or jump_step <= 1:
                    if new_diff < 0:
                        while (CURRENT_FRAME + drift - process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v) < 0:
                            SPIKE_STEP_TARGET -= 1
                        SPIKE_STEP_TARGET += 1
                    else:
                        while (CURRENT_FRAME + drift - process_post_event(SPIKE_STEP_TARGET, mmap_indices, mmap_xy, mmap_data).min_frame_v) >= 0:
                            SPIKE_STEP_TARGET += 1
                    run_flag = False
                else:
                    print(jump_step, sign)

                    if jump_step <= 1:
                        run_flag = False
                    jump_step = jump_step//2
                    if new_diff < 0:
                        sign = -1
                    else:
                        sign = 1
            event_index = SPIKE_STEP_TARGET

    class TimeBar:

        def __init__(self, m_frame, y):
            self.max_frame = m_frame
            self.c_frame = 0
            self.r = 1
            self.g = 1
            self.b = 1
            self.y = 1
            self.r_s = 1
            self.g_s = 1
            self.b_s = 1
            self.alpha = 0.5
            self.d_alpha = 0.05
            self.w = 1
            self.dy = y
            self.w_rem = 0
            self.label = pyglet.text.Label(
                "",
                font_name='Lucida Console',
                font_size=10,
                x=-10000, y=-10000,
                color=(int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), int(self.alpha*255)),
                anchor_x='center',
                anchor_y='bottom')
            self.SELECTING = False
            self.sel_string = ""

        def update(self, c_frame, w_rem):
            self.c_frame = c_frame
            self.w = max(1, int(window.p_width * 0.8))
            self.y = int(self.dy)
            self.w_rem = w_rem
            self.label.color = (int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), 255)#int(self.alpha*255))
            self.label.y = int(1.5*self.y)
            self.label.x = window.p_width//2
            if not self.SELECTING:
                self.label.text = str(c_frame) + '/' + str(self.max_frame)
            else:
                self.label.text = ">" + self.sel_string + '/' + str(self.max_frame)

        def draw(self):
            pyglet.gl.glBegin(pyglet.gl.GL_LINES)
            pyglet.gl.glColor4f(BASE_R, BASE_G, BASE_B, self.alpha)
            pyglet.gl.glVertex2f((window.p_width//2) - (self.w//2), self.y)
            pyglet.gl.glVertex2f((window.p_width//2) + (self.w//2), self.y)
            pyglet.gl.glEnd()

            plx = int((self.c_frame/self.max_frame) * self.w) + ((window.p_width//2) - (self.w//2))

            pyglet.gl.glColor4f(BG_R, BG_G, BG_B, self.alpha)
            circle_og_fill_nozoom(plx, self.y, 10)
            pyglet.gl.glColor4f(BASE_R, BASE_G, BASE_B, self.alpha)
            circle_og_fill_nozoom(plx, self.y, 9)

            pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
            pyglet.gl.glColor3f(BG_R, BG_G, BG_B)

            pyglet.gl.glVertex2f(self.label.x - self.label.content_width//2 - 2, self.label.y - 2)
            pyglet.gl.glVertex2f(self.label.x + self.label.content_width//2 + 2, self.label.y - 2)
            pyglet.gl.glVertex2f(self.label.x + self.label.content_width//2 + 2, self.label.y + self.label.content_height + 2)
            pyglet.gl.glVertex2f(self.label.x - self.label.content_width//2 - 2, self.label.y + self.label.content_height + 2)

            pyglet.gl.glEnd()

            pyglet.gl.glBegin(pyglet.gl.GL_LINES)
            pyglet.gl.glColor3f(BASE_R, BASE_G, BASE_B)

            pyglet.gl.glVertex2f(self.label.x - self.label.content_width//2 - 2, self.label.y - 2)
            pyglet.gl.glVertex2f(self.label.x + self.label.content_width//2 + 2, self.label.y - 2)

            pyglet.gl.glVertex2f(self.label.x + self.label.content_width//2 + 2, self.label.y - 2)
            pyglet.gl.glVertex2f(self.label.x + self.label.content_width//2 + 2, self.label.y + self.label.content_height + 2)

            pyglet.gl.glVertex2f(self.label.x + self.label.content_width//2 + 2, self.label.y + self.label.content_height + 2)
            pyglet.gl.glVertex2f(self.label.x - self.label.content_width//2 - 2, self.label.y + self.label.content_height + 2)

            pyglet.gl.glVertex2f(self.label.x - self.label.content_width//2 - 2, self.label.y + self.label.content_height + 2)
            pyglet.gl.glVertex2f(self.label.x - self.label.content_width//2 - 2, self.label.y - 2)

            pyglet.gl.glEnd()

            self.label.draw()

    class InfoWindow:

        def __init__(self, x, y, w, h, title, titleheight):

            self.fontsize = 12
            self.titlesize = 18

            self.alpha = 0
            self.d_alpha = 1

            self.LOCKED = False
            self.SELECTED = False

            self.title = title
            self.x = x
            self.y = y
            self.w = w + window_radius * 2
            self.h = h + window_radius + max(titleheight, window_radius)
            self.titleheight = max(titleheight, window_radius)
            self.t_label = pyglet.text.Label(
                self.title,
                font_name='Lucida Console',
                font_size=self.titlesize,
                x=(self.x - (self.w/2)), y=(self.y - (self.titleheight/2)),
                color=(int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), int(self.alpha*255)),
                anchor_x='center',
                anchor_y='center')
            self.labels = []
            for i in range(10):
                self.labels.append(pyglet.text.Label(
                    '',
                    font_name='Lucida Console',
                    font_size=self.fontsize,
                    x=-10000, y=-10000,
                    color=(int(BASE_R * 255), int(BASE_G * 255), int(BASE_B * 255), int(self.alpha * 255)),
                    anchor_x='left',
                    anchor_y='top'))
            self.minsize = 200
            self.bdw = 2

            self.event = None
            self.event_info = [0, 0, [0], [], 0]
            self.channel = None
            self.channel_event_fontsize = 10
            self.channel_event_fontsize_2 = 14

            self.channel_border = 4

            self.chn_event_labels = []

            for i in range(STEPS_TO_SEE_MAX):
                self.chn_event_labels.append(pyglet.text.Label("",
                                                          font_name='Lucida Console',
                                                          font_size=self.channel_event_fontsize,
                                                          x=0, y=0,
                                                          anchor_x='left', anchor_y='top'))

            self.mxh = 0
            self.evh = 0

        def select(self, x, y):
            self.SELECTED = True
            num_events = len(self.channel.events) if self.channel is not None else self.evh
            print('sel', num_events)
            for i in range(num_events):
                lbl = self.chn_event_labels[i]
                if self.event is not None:
                    lbl.x += self.x - self.w + self.w//2 - int(self.event_info[2][0])//2
                    lbl.y += self.y - self.h + int(self.event_info[4])
                print(x, y, lbl.x, lbl.content_width, lbl.y, lbl.content_height)
                if lbl.x < x < lbl.x + lbl.content_width and lbl.y - lbl.content_height < y < lbl.y:
                    if self.channel is not None:
                        widget_jump(self.channel.events[i][1].get_x(), self.channel.events[i][1].get_y())
                    elif self.event is not None:
                        chn = channels_linear[int(lbl.text.split("(")[0].strip())]
                        widget_jump(chn.x, chn.y)
                if self.event is not None:
                    lbl.x -= self.x - self.w + self.w//2 - int(self.event_info[2][0])//2
                    lbl.y -= self.y - self.h + int(self.event_info[4])

        def unselect(self):
            self.SELECTED = False

        def text_set(self, m_w, m_w_2, m_w_3):
            max_size = m_w
            if self.channel is not None:
                max_size = max(m_w, m_w_2)
            if self.event is not None:
                max_size = max(m_w, m_w_3)
            for i in range(self.mxh):
                if self.labels[i].content_width > max_size:
                    max_size = self.labels[i].content_width
            max_height = self.mxh * 2 * self.fontsize
            self.w = max(max_size, self.minsize) + 2*self.bdw
            self.h = max_height + self.titleheight + 2*self.bdw
            c_h = self.titleheight
            for i in self.labels:
                i.fontsize = self.fontsize
                i.x = self.x - self.w + self.bdw
                i.y = self.y - c_h - self.bdw
                c_h += 2*self.fontsize

        event_sides = 16

        def get_event_bounds(self):
            if self.event is None:
                return [0, 0, 0, 0]
            min_x = self.event.get_x()
            max_x = self.event.get_x()
            min_y = self.event.get_y()
            max_y = self.event.get_y()
            for e in self.event.data:
                x = channels_linear[e[0]].x
                y = channels_linear[e[0]].y
                if min_x is None:
                    min_x = x
                else:
                    if x < min_x:
                        min_x = x
                if min_y is None:
                    min_y = y
                else:
                    if y < min_y:
                        min_y = y
                if max_x is None:
                    max_x = x
                else:
                    if x > max_x:
                        max_x = x
                if max_y is None:
                    max_y = y
                else:
                    if y > max_y:
                        max_y = y
            bounds = [min_x, min_y, max_x, max_y]
            self.event_info[0] = self.event.get_x() - bounds[0]
            self.event_info[1] = self.event.get_y() - bounds[1]
            self.event_info[2] = [bounds[2] - bounds[0], bounds[3] - bounds[1]]
            event_max_h = window.p_height//2
            event_max_w = window.p_width//2
            scalar = 1
            if self.event_info[2][1] * (event_max_w/self.event_info[2][0]) > event_max_h:
                scalar = event_max_h/self.event_info[2][1]
            else:
                scalar = event_max_w/self.event_info[2][0]
            if scalar > 1:
                scalar = 1
            self.event_info[0] *= scalar
            self.event_info[1] *= scalar
            self.event_info[2][0] *= scalar
            self.event_info[2][1] *= scalar
            self.evh = 0
            for cnt, e in enumerate(self.event.data):
                x = channels_linear[e[0]].x - bounds[0]
                y = channels_linear[e[0]].y - bounds[1]
                x *= scalar
                y *= scalar
                self.event_info[3].append([x, y])
                self.chn_event_labels[cnt].text = str(e[0]) + "(" + str(e[1] - self.event.min_frame_v) + ")"
                self.chn_event_labels[cnt].x = x + 1
                self.chn_event_labels[cnt].y = y + 1
                self.chn_event_labels[cnt].color = (0, 0, 0, 255)
                self.chn_event_labels[cnt].font_size = self.channel_event_fontsize_2
                if self.chn_event_labels[cnt].x + self.chn_event_labels[cnt].content_width + 1 > self.event_info[2][0]:
                    self.event_info[2][0] = self.chn_event_labels[cnt].x + self.chn_event_labels[cnt].content_width + 1
                _ = int(self.event_info[2][1])
                if self.chn_event_labels[cnt].y - self.chn_event_labels[cnt].content_height - 1 < 0:
                    self.event_info[4] = min(self.chn_event_labels[cnt].y - self.chn_event_labels[cnt].content_height - 1, self.event_info[4])
                    self.event_info[2][1] = max(self.event_info[2][1], _ - (self.chn_event_labels[cnt].y - self.chn_event_labels[cnt].content_height - 1))
                self.evh += 1
            self.event_info[4] = -self.event_info[4]

        def adjust(self, dx, dy):
            self.x += dx
            self.y += dy

        def clear(self):
            self.channel = None
            self.event = None
            self.title = ""
            for i in self.labels:
                i.x = -10000
                i.y = -10000
                i.text = ''
            self.mxh = 0
            self.event_info = [0, 0, [0], [], 0]

        def set_channel(self, channel, counts):
            self.clear()
            self.channel = channel
            self.title = "Channel " + str(channel.index)
            self.labels[0].text = 'x: ' + str(channel.x)
            self.labels[1].text = 'y: ' + str(channel.y)
            self.labels[2].text = 'Spike freq /min: ' + str(np.round(counts[self.channel.index]/N_FRAMES*60,2))
            self.labels[3].text = 'Spike count: ' + str(counts[self.channel.index])
            for label in self.chn_event_labels:
                label.font_size = self.channel_event_fontsize
            self.mxh = 4

        def set_event(self, event):
            self.clear()
            self.event = event
            self.title = "Event " + str(hex(event.identification))
            self.labels[0].text = 'channels: ' + str([s[0] for s in event.data])[1:][:-1]
            self.labels[1].text = 'min. frame: ' + str(self.event.min_frame_v) + " (" + str(self.event.max_frame_v - self.event.min_frame_v) + ")"
            self.mxh = 2
            self.get_event_bounds()

        def draw(self):

            if self.event is None and self.channel is None:
                self.alpha -= self.d_alpha
                self.alpha = max(0, min(1, self.alpha))
            else:
                self.alpha += self.d_alpha
                self.alpha = max(0, min(1, self.alpha))

            if self.alpha == 0:
                return

            self.t_label.font_size = self.titlesize
            self.t_label.x = (self.x - (self.w/2))
            self.t_label.y = (self.y - (self.titleheight/2))
            self.t_label.text = self.title
            self.t_label.color = (int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), int(self.alpha*255))

            num_events = 0
            chn_indices = []

            self.text_set(self.t_label.content_width, STEPS_TO_SEE * 2, int(self.event_info[2][0]))

            if self.event is not None:
                num_events = self.evh
                self.h += int(self.event_info[2][1])
                event_width = int(self.event_info[2][0])

            if self.channel is not None and STEPS_TO_SEE > 1:
                self.channel.relight()
                chn_queue = []
                # For clarity
                val = -AF_DATA[(0 + mod_count) % STEPS_TO_SEE_MAX][self.channel.index] + self.y - self.h
                min_height = None
                max_height = None
                if STEPS_TO_SEE == STEPS_TO_SEE_MAX:
                    min_height = val
                    max_height = val
                chn_queue.append(val)
                for i in range(1, STEPS_TO_SEE_MAX, 1):
                    val = -AF_DATA[(i + mod_count) % STEPS_TO_SEE_MAX][self.channel.index] + self.y - self.h
                    if min_height is None and i >= STEPS_TO_SEE_MAX - STEPS_TO_SEE:
                        min_height = val
                        max_height = val
                    elif min_height is not None:
                        if val > max_height:
                            max_height = val
                        if val < min_height:
                            min_height = val
                    chn_queue.append(val)
                diff = abs(min_height - max_height)
                for tlt in self.channel.events:
                    if STEPS_TO_SEE_MAX - STEPS_TO_SEE <= STEPS_TO_SEE_MAX - CURRENT_FRAME + tlt[0] < STEPS_TO_SEE_MAX:
                        b_index = STEPS_TO_SEE - CURRENT_FRAME + tlt[0]
                        self.chn_event_labels[num_events].text = str(hex(tlt[1].identification))
                        self.chn_event_labels[num_events].font_size = self.channel_event_fontsize
                        self.chn_event_labels[num_events].x = int(b_index * 2 + (self.x - self.w / 2) - (STEPS_TO_SEE / 2) * 2)
                        self.chn_event_labels[num_events].y = self.y - self.h
                        self.chn_event_labels[num_events].color = (int(tlt[1].r * 255), int(tlt[1].g * 255), int(tlt[1].b * 255), int(self.alpha * 255))
                        chn_indices.append(STEPS_TO_SEE_MAX - CURRENT_FRAME + tlt[0])
                        num_events += 1
                self.h += diff + (self.channel_border * 2) + (2 * self.channel_event_fontsize)
                for i in range(len(chn_queue)):
                    chn_queue[i] -= max_height - self.y + self.h - diff - self.channel_border
                w_queue = deque()
                prv = int((self.x - self.w / 2) - (STEPS_TO_SEE / 2) * 2)
                for j in range(STEPS_TO_SEE):
                    w_queue.append(prv)
                    prv += 2
                for j in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE):
                    w_queue.appendleft(-1)

            pyglet.gl.glColor4f(BG_R, BG_G, BG_B, self.alpha)

            circle_quarter_filled(self.x - window_radius, self.y - window_radius, window_radius, 0)
            circle_quarter_filled(self.x - self.w + window_radius, self.y - window_radius, window_radius, 1)
            circle_quarter_filled(self.x - window_radius, self.y - self.h + window_radius, window_radius, 3)
            circle_quarter_filled(self.x - self.w + window_radius, self.y - self.h + window_radius, window_radius, 2)

            pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)

            pyglet.gl.glVertex2f(self.x - self.w + window_radius, self.y)
            pyglet.gl.glVertex2f(self.x - window_radius, self.y)
            pyglet.gl.glVertex2f(self.x - window_radius, self.y - window_radius)
            pyglet.gl.glVertex2f(self.x - self.w + window_radius, self.y - window_radius)

            pyglet.gl.glEnd()

            pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)

            pyglet.gl.glVertex2f(self.x - self.w, self.y - window_radius)
            pyglet.gl.glVertex2f(self.x, self.y - window_radius)
            pyglet.gl.glVertex2f(self.x, self.y - self.titleheight)
            pyglet.gl.glVertex2f(self.x - self.w, self.y - self.titleheight)

            pyglet.gl.glEnd()

            pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)

            pyglet.gl.glVertex2f(self.x - self.w, self.y - self.titleheight)
            pyglet.gl.glVertex2f(self.x, self.y - self.titleheight)
            pyglet.gl.glVertex2f(self.x, self.y - self.h + window_radius)
            pyglet.gl.glVertex2f(self.x - self.w, self.y - self.h + window_radius)

            pyglet.gl.glEnd()

            pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)

            pyglet.gl.glVertex2f(self.x - self.w + window_radius, self.y - self.h + window_radius)
            pyglet.gl.glVertex2f(self.x - window_radius, self.y - self.h + window_radius)
            pyglet.gl.glVertex2f(self.x - window_radius, self.y - self.h)
            pyglet.gl.glVertex2f(self.x - self.w + window_radius, self.y - self.h)

            pyglet.gl.glEnd()

            pyglet.gl.glColor4f(BASE_R, BASE_G, BASE_B, self.alpha)

            circle_quarter(self.x - window_radius, self.y - window_radius, window_radius, 0)
            circle_quarter(self.x - self.w + window_radius, self.y - window_radius, window_radius, 1)
            circle_quarter(self.x - window_radius, self.y - self.h + window_radius, window_radius, 3)
            circle_quarter(self.x - self.w + window_radius, self.y - self.h + window_radius, window_radius, 2)

            pyglet.gl.glBegin(pyglet.gl.GL_LINES)
            #pyglet.gl.glColor4f(BASE_R, BASE_G, BASE_B, 1)

            pyglet.gl.glVertex2f(self.x - self.w + window_radius, self.y)
            pyglet.gl.glVertex2f(self.x - window_radius, self.y)

            pyglet.gl.glVertex2f(self.x - self.w + window_radius, self.y - self.h)
            pyglet.gl.glVertex2f(self.x - window_radius, self.y - self.h)

            pyglet.gl.glVertex2f(self.x - self.w, self.y - self.h + window_radius)
            pyglet.gl.glVertex2f(self.x - self.w, self.y - window_radius)

            pyglet.gl.glVertex2f(self.x, self.y - self.h + window_radius)
            pyglet.gl.glVertex2f(self.x, self.y - window_radius)

            pyglet.gl.glVertex2f(self.x - self.w, self.y - self.titleheight)
            pyglet.gl.glVertex2f(self.x, self.y - self.titleheight)

            pyglet.gl.glEnd()

            self.t_label.draw()

            if self.event is not None:
                t_r, t_g, t_b = self.event.r, self.event.g, self.event.b
                pyglet.gl.glColor4f(t_r, t_g, t_b, self.alpha)
                for point in self.event_info[3]:
                    pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                    pyglet.gl.glVertex2f(self.x - self.w + self.event_info[0] + self.w//2 - event_width//2, self.y - self.h + self.event_info[1] + self.event_info[4])
                    pyglet.gl.glVertex2f(self.x - self.w + point[0] + self.w//2 - event_width//2, self.y - self.h + point[1] + self.event_info[4])
                    pyglet.gl.glEnd()
                circle_og_fill_nozoom(self.x - self.w + self.event_info[0] + self.w//2 - event_width//2, self.y - self.h + self.event_info[1] + self.event_info[4], self.event_sides//2)
                pyglet.gl.glColor4f(BG_R, BG_G, BG_B, self.alpha)
                circle_og_fill_nozoom(self.x - self.w + self.event_info[0] + self.w//2 - event_width//2, self.y - self.h + self.event_info[1] + self.event_info[4], self.event_sides//3)

            if self.channel is not None and STEPS_TO_SEE > 1:
                window.draw_multiline_2(w_queue, chn_queue, self.channel.l_queue)

            for i in self.labels:
                i.color = (int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), int(self.alpha*255))
                i.draw()

            if self.channel is not None:
                for tlt in range(num_events):

                    t_r, t_g, t_b = self.chn_event_labels[tlt].color[0]/255, self.chn_event_labels[tlt].color[1]/255, self.chn_event_labels[tlt].color[2]/255

                    pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                    pyglet.gl.glColor4f(t_r, t_g, t_b, self.alpha)

                    pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x, self.chn_event_labels[tlt].y - self.chn_event_labels[tlt].content_height)
                    pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x, chn_queue[chn_indices[tlt]] + 5)

                    pyglet.gl.glEnd()

                    if self.chn_event_labels[tlt].x + self.chn_event_labels[tlt].content_width < self.x:
                        pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
                        pyglet.gl.glColor4f(t_r, t_g, t_b, self.alpha)

                        pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x - 1, self.chn_event_labels[tlt].y + 1)
                        pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x + self.chn_event_labels[tlt].content_width + 1,
                                             self.chn_event_labels[tlt].y + 1)
                        pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x + self.chn_event_labels[tlt].content_width + 1,
                                             self.chn_event_labels[tlt].y - 1 - self.chn_event_labels[tlt].content_height)
                        pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x - 1,
                                             self.chn_event_labels[tlt].y - self.chn_event_labels[tlt].content_height - 1)

                        pyglet.gl.glEnd()
                        self.chn_event_labels[tlt].color = (0, 0, 0, int(self.alpha * 255))
                        self.chn_event_labels[tlt].draw()

            elif self.event is not None:
                t_r, t_g, t_b = self.event.r, self.event.g, self.event.b

                for tlt in range(num_events):
                    self.chn_event_labels[tlt].x = self.chn_event_labels[tlt].x + self.x - self.w + self.w//2 - event_width//2
                    self.chn_event_labels[tlt].y = self.chn_event_labels[tlt].y + self.y - self.h + int(self.event_info[4])
                    #self.chn_event_labels[tlt].color = (int(t_r * 255), int(t_g * 255), int(t_b * 255), int(self.alpha * 255))

                    pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
                    pyglet.gl.glColor4f(t_r, t_g, t_b, self.alpha)

                    pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x - 1, self.chn_event_labels[tlt].y + 1)
                    pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x + self.chn_event_labels[tlt].content_width + 1,
                                         self.chn_event_labels[tlt].y + 1)
                    pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x + self.chn_event_labels[tlt].content_width + 1,
                                         self.chn_event_labels[tlt].y - 1 - self.chn_event_labels[tlt].content_height)
                    pyglet.gl.glVertex2f(self.chn_event_labels[tlt].x - 1,
                                         self.chn_event_labels[tlt].y - self.chn_event_labels[tlt].content_height - 1)

                    pyglet.gl.glEnd()

                    self.chn_event_labels[tlt].draw()
                    self.chn_event_labels[tlt].x -= self.x - self.w + self.w//2 - event_width//2
                    self.chn_event_labels[tlt].y -= self.y - self.h + + int(self.event_info[4])

    class VCanvas(pyglet.window.Window):

        def __init__(self):
            config = pyglet.gl.Config(double_buffer=False)
            super().__init__(config=config, fullscreen=False, resizable=True)
            self.batch = pyglet.graphics.Batch()
            self.p_width = self.width
            self.p_height = self.height

        def draw_line(self, x1, y1, x2, y2):
            self.batch.add(2, pyglet.gl.GL_LINES, None, ("v2f", (x1, y1, x2, y2)))

        def scroll_x(self, x_sc):
            global ox
            ox += x_sc

        def scroll_y(self, y_sc):
            global oy
            oy -= y_sc

        def draw_multiline_2(self, xs, ys, mem):

            if render_mode == 0:
                pass

            elif render_mode == 1:
                pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                if THEME == 0:
                    for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                        pyglet.gl.glColor3f(BASE_R, BASE_G-mem[i]/100, BASE_B-mem[i]/100)
                        pyglet.gl.glVertex2f(xs[i], ys[i])
                        pyglet.gl.glVertex2f(xs[i + 1], ys[i + 1])
                else:
                    for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                        pyglet.gl.glColor3f(BASE_R+mem[i]/100, BASE_G, BASE_B)
                        pyglet.gl.glVertex2f(xs[i], ys[i])
                        pyglet.gl.glVertex2f(xs[i + 1], ys[i + 1])
                pyglet.gl.glEnd()

        def draw_multiline_2_decay(self, xs, ys, mem):

            if render_mode == 0:
                pass

            elif render_mode == 1:
                pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                if THEME == 0:
                    for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                        pyglet.gl.glColor3f(BASE_R, BASE_G-mem[i]/100, BASE_B-mem[i]/100)
                        if i - STEPS_TO_SEE_MAX + STEPS_TO_SEE <= STEPS_TO_SEE//4:
                            brt = (i - STEPS_TO_SEE_MAX + STEPS_TO_SEE) / max((STEPS_TO_SEE // 4), 1)
                            pyglet.gl.glColor4f(BASE_R, BASE_G - mem[i] / 100, BASE_B - mem[i] / 100, brt)
                        pyglet.gl.glVertex2f(xs[i], ys[i])
                        pyglet.gl.glVertex2f(xs[i + 1], ys[i + 1])
                else:
                    for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                        pyglet.gl.glColor3f(BASE_R+mem[i]/100, BASE_G, BASE_B)
                        if i - STEPS_TO_SEE_MAX + STEPS_TO_SEE <= STEPS_TO_SEE//4:
                            brt = (i - STEPS_TO_SEE_MAX + STEPS_TO_SEE) / max((STEPS_TO_SEE // 4), 1)
                            pyglet.gl.glColor4f(BASE_R+mem[i]/100, BASE_G, BASE_B, brt)
                        pyglet.gl.glVertex2f(xs[i], ys[i])
                        pyglet.gl.glVertex2f(xs[i + 1], ys[i + 1])
                pyglet.gl.glEnd()

        def draw_multiline(self, xs, ys):

            if render_mode == 0:
                pyglet.gl.glColor3f(BASE_R, BASE_G, BASE_B)
                for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                    self.draw_line(xs[i], ys[i], xs[i+1], ys[i+1])

            elif render_mode == 1:
                pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                pyglet.gl.glColor3f(BASE_R, BASE_G, BASE_B)
                for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                    pyglet.gl.glVertex2f(xs[i], ys[i])
                    pyglet.gl.glVertex2f(xs[i + 1], ys[i + 1])
                pyglet.gl.glEnd()

        def draw_multiline_decay(self, xs, ys):

            if render_mode == 0:
                pyglet.gl.glColor3f(BASE_R, BASE_G, BASE_B)
                for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                    self.draw_line(xs[i], ys[i], xs[i+1], ys[i+1])

            elif render_mode == 1:
                pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                for i in range(STEPS_TO_SEE_MAX - STEPS_TO_SEE, STEPS_TO_SEE_MAX - 1, 1):
                    pyglet.gl.glColor3f(BASE_R, BASE_G, BASE_B)
                    if i - STEPS_TO_SEE_MAX + STEPS_TO_SEE <= STEPS_TO_SEE // 4:
                        brt = (i - STEPS_TO_SEE_MAX + STEPS_TO_SEE) / max((STEPS_TO_SEE // 4), 1)
                        pyglet.gl.glColor4f(BASE_R, BASE_G, BASE_B, brt)
                    pyglet.gl.glVertex2f(xs[i], ys[i])
                    pyglet.gl.glVertex2f(xs[i + 1], ys[i + 1])
                pyglet.gl.glEnd()

        def batch_draw(self):
            self.batch.draw()
            self.batch = pyglet.graphics.Batch()

    window = VCanvas()

    half_w = int(window.p_width/2)
    half_h = int(window.p_height/2)

    platform = pyglet.window.get_platform()
    display = platform.get_default_display()
    pyglet.gl.glEnable(pyglet.gl.GL_BLEND)
    pyglet.gl.glBlendFunc(pyglet.gl.GL_SRC_ALPHA, pyglet.gl.GL_ONE_MINUS_SRC_ALPHA)

    folder_name = folder_name.strip()

    if not os.path.isdir(os.path.join("gen", folder_name)):
        os.mkdir(os.path.join("gen", folder_name))

    root.destroy()

    try:
        path_sizes = os.path.join('gen', folder_name, 'meta_sizes')
        path_indices = os.path.join('gen', folder_name, 'meta_indices')
        path_data = os.path.join('gen', folder_name, 'meta_data')
        path_xy = os.path.join('gen', folder_name, 'meta_xypos')
        path_totals = os.path.join('gen', folder_name, 'meta_totals')
        path_sc = os.path.join('gen', folder_name, 'meta_spikecount')
        path_averages = os.path.join('gen', folder_name, 'meta_averages')
    except Exception as e:
        errnotify.showerror('OS Path Joining Exception', 'Something went wrong during OS path conversion: ' + str(e))
        exit()

    to_pass = [path_sizes, path_indices, path_data, path_xy, path_totals, path_sc, path_averages]

    to_load = True
    if skip_loading:
        to_load = False

    prev_ld = 0
    phase = 0
    dots = 0

    label = pyglet.text.Label("",
                              font_name='Lucida Console',
                              font_size=20,
                              x=10, y=10,
                              color=(int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), 255),
                              anchor_x='left', anchor_y='bottom')

    if to_load:
        load_process = mul.Process(target=lr.textfile_processor, args=(load_end_flag, shape1, shape2, shape3, file_proc, squares_processed.shape, is_biocam, len(probe.positions), avg, to_pass))
        load_process.start()

    while to_load and load_end_flag.value < 100.5:
        window.dispatch_events()
        window.clear()
        if load_end_flag.value < 0:
            dots += 1
            dots %= 60
            dd = ""
            for i in range(int(dots/20)):
                dd += "."
            label.text = "Counting lines." + dd
        else:
            a_str = ", querying indices" if phase == 1 else (
            ", modifying indices" if phase == 2 else ", writing indices")
            label.text = "Loading" + a_str + " (" + str(phase) + "/3): " + str(round(load_end_flag.value, 2)) + "%"

        if prev_ld > load_end_flag.value:
            phase += 1
        prev_ld = load_end_flag.value
        label.draw()
        pyglet.gl.glFlush()
        time.sleep(0.01)

    if to_load:
        window.dispatch_events()
        window.clear()
        label.text = "Loading: 100-ish%. Almost done!"
        prev_ld = load_end_flag.value
        label.draw()
        pyglet.gl.glFlush()

        load_process.terminate()

    try:
        bool_sum = True
        bool_sum = bool_sum and os.path.exists(path_sizes)
        bool_sum = bool_sum and os.path.exists(path_indices)
        bool_sum = bool_sum and os.path.exists(path_data)
        bool_sum = bool_sum and os.path.exists(path_xy)
        bool_sum = bool_sum and os.path.exists(path_totals)
        bool_sum = bool_sum and os.path.exists(path_sc)
        if not bool_sum:
            errnotify.showerror('Missing Generated File Exception', 'One of the required files from gen/' + folder_name + ' is missing. Please re-run the generation process.')
            exit()
    except Exception:
        errnotify.showerror('OS File Checking Exception', 'Something went wrong during OS file checking: ' + str(e))
        exit()

    mmap_sizes = np.memmap(path_sizes, mode='r', shape=(3,), dtype='int64')

    shape1.value = mmap_sizes[0]
    shape2.value = mmap_sizes[1]
    shape3.value = mmap_sizes[2]

    mmap_indices = np.memmap(path_indices, mode='r', shape=(2, shape1.value), dtype='int32')
    mmap_data = np.memmap(path_data, mode='r', shape=(3, shape2.value), dtype='int32')
    mmap_xy = np.memmap(path_xy, mode='r', shape=(3, shape3.value), dtype='float32')
    mmap_totals = np.memmap(path_totals, mode='r', shape=squares_processed.shape, dtype='float32')
    mmap_totals_proc = gen_smoothed_arrays(mmap_totals)
    mmap_spikecount = np.memmap(path_sc, mode='r', shape=len(probe.positions), dtype='int64')

    post_processed = True

    CURRENT_FRAME = 0

    # 0: Batch rendering,  1: glVertex rendering,  2: indexed list rendering
    render_mode = 1

    ox = 0
    oy = 0

    dx = 0
    dy = 0
    dzoom = 0

    zoom = 1.0

    MEMORY_MODEL = 0  # 0: memmap, 1: memmap + main mem.

    MAIN_MEMORY_WINDOW = 50000
    MAIN_MEMORY_DATA = probe.Read(0, MAIN_MEMORY_WINDOW).reshape((MAIN_MEMORY_WINDOW, data_base.shape[0]))
    MAIN_MEMORY_COUNT = 0
    MAIN_MEMORY_COUNT_PLUS = 0

    AF_DATA = []
    for i in range(STEPS_TO_SEE_MAX):
        AF_DATA.append(np.zeros(data_base.shape))

    mod_count = 0

    def circle_quarter(x, y, radius, quarter):
        iterations = max(int(radius * math.pi), 1)
        s = math.sin(2 * math.pi / iterations)
        c = math.cos(2 * math.pi / iterations)

        dx, dy = radius, 0
        if quarter == 1:
            dx, dy = 0, radius
        elif quarter == 2:
            dx, dy = -radius, 0
        elif quarter == 3:
            dx, dy = 0, -radius

        pyglet.gl.glBegin(pyglet.gl.GL_LINES)
        for i in range(iterations//4 + 1):
            pyglet.gl.glVertex2f(x + dx, y + dy)
            dx, dy = (dx * c - dy * s), (dy * c + dx * s)
            pyglet.gl.glVertex2f(x + dx, y + dy)
        pyglet.gl.glEnd()

    def circle_quarter_filled(x, y, radius, quarter):
        iterations = max(int(radius * math.pi), 1)
        s = math.sin(2 * math.pi / iterations)
        c = math.cos(2 * math.pi / iterations)

        dx, dy = radius, 0
        if quarter == 1:
            dx, dy = 0, radius
        elif quarter == 2:
            dx, dy = -radius, 0
        elif quarter == 3:
            dx, dy = 0, -radius

        pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
        pyglet.gl.glVertex2f(x, y)
        for i in range(iterations//4 + 1):
            pyglet.gl.glVertex2f(x + dx, y + dy)
            dx, dy = (dx * c - dy * s), (dy * c + dx * s)
        pyglet.gl.glEnd()

    def circle(x, y, radius):
        radius = radius * zoom
        iterations = max(int(radius * math.pi), 1)
        s = math.sin(2 * math.pi / iterations)
        c = math.cos(2 * math.pi / iterations)

        dx, dy = radius, 0

        s_x = (x + ox) * zoom + half_w
        s_y = (y + oy) * zoom + half_h

        pyglet.gl.glBegin(pyglet.gl.GL_LINES)
        for i in range(iterations + 1):
            pyglet.gl.glVertex2f(s_x + dx, s_y + dy)
            dx, dy = (dx * c - dy * s), (dy * c + dx * s)
            pyglet.gl.glVertex2f(s_x + dx, s_y + dy)
        pyglet.gl.glEnd()

    def circle_fill(x, y, radius):
        radius = radius * zoom
        iterations = max(int(radius * math.pi), 1)
        s = math.sin(2 * math.pi / iterations)
        c = math.cos(2 * math.pi / iterations)

        dx, dy = radius, 0

        s_x = (x + ox) * zoom + half_w
        s_y = (y + oy) * zoom + half_h

        pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
        pyglet.gl.glVertex2f(s_x, s_y)
        for i in range(iterations + 1):
            pyglet.gl.glVertex2f(s_x + dx, s_y + dy)
            dx, dy = (dx * c - dy * s), (dy * c + dx * s)
            pyglet.gl.glVertex2f(s_x + dx, s_y + dy)
        pyglet.gl.glEnd()

    def circle_og(x, y, radius):
        radius = radius * zoom
        iterations = max(int(radius * math.pi), 1)
        s = math.sin(2 * math.pi / iterations)
        c = math.cos(2 * math.pi / iterations)

        dx, dy = radius, 0

        pyglet.gl.glBegin(pyglet.gl.GL_LINES)
        for i in range(iterations + 1):
            pyglet.gl.glVertex2f(x + dx, y + dy)
            dx, dy = (dx * c - dy * s), (dy * c + dx * s)
            pyglet.gl.glVertex2f(x + dx, y + dy)
        pyglet.gl.glEnd()

    def circle_og_fill(x, y, radius):
        radius = radius * zoom
        iterations = max(int(radius * math.pi), 1)
        s = math.sin(2 * math.pi / iterations)
        c = math.cos(2 * math.pi / iterations)

        dx, dy = radius, 0

        pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
        pyglet.gl.glVertex2f(x, y)
        for i in range(iterations + 1):
            pyglet.gl.glVertex2f(x + dx, y + dy)
            dx, dy = (dx * c - dy * s), (dy * c + dx * s)
            pyglet.gl.glVertex2f(x + dx, y + dy)
        pyglet.gl.glEnd()

    def circle_og_fill_nozoom(x, y, radius):
        iterations = max(int(radius * math.pi), 1)
        s = math.sin(2 * math.pi / iterations)
        c = math.cos(2 * math.pi / iterations)

        dx, dy = radius, 0

        pyglet.gl.glBegin(pyglet.gl.GL_TRIANGLE_FAN)
        pyglet.gl.glVertex2f(x, y)
        for i in range(iterations + 1):
            pyglet.gl.glVertex2f(x + dx, y + dy)
            dx, dy = (dx * c - dy * s), (dy * c + dx * s)
            pyglet.gl.glVertex2f(x + dx, y + dy)
        pyglet.gl.glEnd()

    def add_channels(lst):
        for i, tup in enumerate(lst):
            if is_biocam:
                h = Channel(i, tup[0] + 100, tup[1] + 100)
            else:
                h = Channel(i, tup[0], tup[1])
            h.append_this(channels, B_WIDTH, B_HEIGHT)
            channels_linear.append(h)
        return channels

    light_decay = 5

    class Channel:

        def __init__(self, index, x, y):
            self.vertex_list = pyglet.graphics.vertex_list(100, 'v2f/stream')
            self.queue = deque() #y-values
            self.l_queue = deque() #brightness-values
            self.x_queue = deque() #x-values
            for i in range(STEPS_TO_SEE_MAX):
                self.x_queue.append(0)
            for i in range(STEPS_TO_SEE_MAX):
                self.l_queue.append(0)
            for i in range(STEPS_TO_SEE_MAX):
                self.queue.append(y + half_h)
            self.bts = 0
            self.l_total = 0
            self.index = index
            self.x = x
            self.y = y
            self.step_xdiv = 2
            self.adjust()
            self.LOADED = False
            self.p_zoom = zoom
            self.p_ox = ox
            self.p_oy = oy
            self.p_frame = 0
            self.to_light = []
            self.events = []

        def trigger_spike(self, frame):
            self.to_light.append(frame)
            self.to_light = [s for s in self.to_light if s >= CURRENT_FRAME - STEPS_TO_SEE_MAX]

        def relight(self):
            self.to_light = [s for s in self.to_light if s >= CURRENT_FRAME - STEPS_TO_SEE_MAX]
            temp_light = []
            for item in self.to_light:
                temp_light.append(item - CURRENT_FRAME + STEPS_TO_SEE_MAX)
            self.l_total = 0
            temp_bts = 0
            for i in range(STEPS_TO_SEE_MAX):
                if i in temp_light:
                    temp_bts = 100
                self.l_total += temp_bts
                self.l_queue[i] = temp_bts
                temp_bts -= light_decay
                temp_bts = max(0, temp_bts)
            self.bts = temp_bts

        def add_event(self, frame, event):
            self.events.append([frame, event])

        def update(self, value):
            for item in self.to_light:
                if item == CURRENT_FRAME:
                    self.bts = 100
            self.to_light = [s for s in self.to_light if s >= CURRENT_FRAME - STEPS_TO_SEE_MAX]
            self.events = [s for s in self.events if s[0] >= CURRENT_FRAME - STEPS_TO_SEE_MAX]
            self.l_queue.append(self.bts)
            self.l_total += self.bts
            self.bts -= light_decay
            self.bts = max(0, self.bts)
            self.l_total -= self.l_queue.popleft()
            if not self.p_zoom == zoom:
                self.up_zoom()
            elif not self.p_oy == oy:
                self.regen(AF_DATA)
            else:
                self.queue.append(half_h + (-value + self.y + oy) * zoom)
                self.queue.popleft()
            if not self.p_ox == ox:
                self.adjust()
            self.untag()

        def untag(self):
            for item in self.events:
                item[1].TAGGED = False

        def recalc(self):
            if not self.p_zoom == zoom:
                self.up_zoom()
            if not self.p_ox == ox:
                self.adjust()
            if not self.p_oy == oy:
                self.regen(AF_DATA)
            self.untag()

        def up_zoom(self):
            self.adjust()
            self.regen(AF_DATA)
            self.p_zoom = zoom

        def regen(self, lst):
            for i in range(0, STEPS_TO_SEE_MAX, 1):
                self.queue[i] = half_h + (-lst[(i + mod_count) % STEPS_TO_SEE_MAX][self.index] + self.y + oy) * zoom
            self.p_oy = oy

        def load(self):
            self.LOADED = True
            self.events = [s for s in self.events if s[0] >= CURRENT_FRAME - STEPS_TO_SEE_MAX]
            self.up_zoom()
            self.relight()

        def unload(self):
            self.LOADED = False

        def adjust(self):
            prv = half_w + int((self.x + ox - (self.step_xdiv * ((STEPS_TO_SEE / 2) + ((STEPS_TO_SEE_MAX - STEPS_TO_SEE))))) * zoom)
            m_prv = 0
            for i, xi in enumerate(self.x_queue):
                self.x_queue[i] = prv + (m_prv * zoom)
                m_prv += self.step_xdiv
            self.p_ox = ox

        def y_adjust(self):
            diff = self.p_oy - oy
            for item in range(len(self.queue)):
                self.queue[item] -= (diff*zoom)
            self.p_oy = oy

        def zoom_adjust(self):
            self.adjust()
            for item in range(len(self.queue)):
                self.queue[item] -= half_h
                self.queue[item] /= self.p_zoom
                self.queue[item] *= zoom
                self.queue[item] += half_h
            self.y_adjust()
            self.p_oy = oy
            self.p_zoom = zoom

        def append_this(self, lst, band_width, band_height):
            try:
                lst[self.x // band_width][self.y // band_height].append(self)
            except IndexError:
                print("SpikeSeer (ERROR): ", self.x // band_width, self.y // band_height, " !!!")

    class Event:

        def __init__(self, x, y, channel_frame, to_process, index, main_channel, IDEN):

            self.mc = int(main_channel)
            self.x = x
            self.y = y
            self.data = channel_frame
            self.r = 1
            self.g = 1
            self.b = 1
            if to_process:
                self.ident = 0
                self.min_frame_v = self.min_frame()
                self.max_frame_v = self.max_frame()
                self.index = index
            self.TAGGED = False
            self.identification = IDEN

        def __str__(self):
            return str(self.x) + ", " + str(self.y) + ", " + str(self.data) + ", " + str(self.min_frame_v) + ", " + str(self.max_frame_v)

        def min_frame(self):
            num = N_FRAMES
            for step in self.data:
                if step[1] < num:
                    num = step[1]
            return num

        def max_frame(self):
            num = 0
            for step in self.data:
                if step[1] > num:
                    num = step[1]
            return num

        def get_x(self):
            if is_biocam:
                return (self.x * GLOBAL_SCALAR) + 100
            return self.x * GLOBAL_SCALAR

        def get_y(self):
            if is_biocam:
                return (self.y * GLOBAL_SCALAR) + 100
            return self.y * GLOBAL_SCALAR

    def process_post_event(index, idfile, xyfile, datafile):

        x = xyfile[0][index]
        y = xyfile[1][index]

        main_channel = xyfile[2][index]
        begin = idfile[0][index]
        end = idfile[1][index]

        data = []

        for num in range(begin, end):
            data.append([datafile[0][num], datafile[1][num], datafile[2][num]])

        data = np.asarray(data)

        return Event(x, y, data, True, index, main_channel, event_index)

    add_channels(MOD_POS)

    time_counter = 0

    paused = False

    class SmartDict(dict):
        def __missing__(self, key):
            return False

    pressed = SmartDict()

    def handle_zoom(scalar):
        global zoom, dzoom
        zoom *= scalar
        if zoom > 2:
            zoom = 2
            dzoom = 0

    info = InfoWindow(window.p_width - 10, window.p_height - 10, 200, 300, "Info", 40)
    info_list = []
    timebar = TimeBar(N_FRAMES, 20)

    def distance(x1, y1, x2, y2):
        diff_x = x1 - x2
        diff_y = y1 - y2
        return math.sqrt((diff_x * diff_x) + (diff_y * diff_y))

    def find_nearest_object(x, y, rad, passive, snd_win=None):

        if info.LOCKED and passive:
            return

        rad /= zoom
        rad = int(rad)
        s_x = int(((x - half_w)/zoom) - ox)
        s_y = int(((y - half_h)/zoom) - oy)
        lx = (s_x - rad)//B_WIDTH
        rx = (s_x + rad)//B_WIDTH
        dy = (s_y - rad)//B_HEIGHT
        uy = (s_y + rad)//B_HEIGHT
        lx = max(0, min(WIDTH + 1, lx))
        rx = max(0, min(WIDTH + 1, rx))
        dy = max(0, min(HEIGHT + 1, dy))
        uy = max(0, min(HEIGHT + 1, uy))

        dist_min = 2 * rad
        t_obj = None
        o_type = -1

        for i in range(lx, rx + 1, 1):
            for j in range(dy, uy + 1, 1):
                for chn in channels[i][j]:
                    chn_x = chn.x
                    chn_y = chn.y
                    dst = distance(s_x, s_y, chn_x, chn_y)
                    if dst < dist_min:
                        dist_min = dst
                        t_obj = chn
                        o_type = 0
                    for event in chn.events:
                        ev = event[1]
                        ev_x = int(ev.x * GLOBAL_SCALAR)
                        ev_y = int(ev.y * GLOBAL_SCALAR)
                        dst = distance(s_x, s_y, ev_x, ev_y)
                        if dst < dist_min:
                            dist_min = dst
                            t_obj = ev
                            o_type = 1

        if snd_win is not None:
            if o_type != -1:
                if o_type == 0:
                    snd_win.set_channel(t_obj, mmap_spikecount)
                elif o_type == 1:
                    snd_win.set_event(t_obj)
            else:
                snd_win.clear()
            return

        if o_type != -1:
            if o_type == 0:
                info.set_channel(t_obj, mmap_spikecount)
            elif o_type == 1:
                info.set_event(t_obj)
        else:
            info.clear()

    @window.event
    def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
        global ox, oy
        if SETTINGS_MODE:
            return
        if buttons & pyglet.window.mouse.LEFT:
            for inf in reversed(info_list):
                if inf.t_label.x - inf.t_label.content_width/2 < x - dx < inf.t_label.x + inf.t_label.content_width/2 and inf.t_label.y - inf.t_label.content_height/2 < y - dy < inf.t_label.y + inf.t_label.content_height/2:
                    inf.x += dx
                    inf.y += dy
                    return
            print(dx, dy, 'm_drag')
            ox += int(dx/zoom)
            oy += int(dy/zoom)

    @window.event
    def on_mouse_motion(x, y, dx, dy):
        if SETTINGS_MODE:
            return
        focus = None

        if focus is None:
            if not info.LOCKED:
                find_nearest_object(x, y, 100, True)

    def full_unload():
        pv_x = int((-ox - (half_w/zoom)) // B_WIDTH) - 1
        pv_x_2 = int(pv_x + ((2 * half_w/zoom) // B_WIDTH)) + 1
        pv_x = max(0, min(WIDTH + 1, pv_x))
        pv_x_2 = max(0, min(WIDTH + 1, pv_x_2))

        pv_y = int((-oy - (half_h/zoom)) // B_HEIGHT)
        pv_y_2 = int(pv_y + ((2 * half_h/zoom) // B_HEIGHT)) + 2
        pv_y = max(0, min(HEIGHT + 1, pv_y))
        pv_y_2 = max(0, min(HEIGHT + 1, pv_y_2))

        for i in range(pv_x, pv_x_2 + 1, 1):
            for j in range(pv_y, pv_y_2 + 1, 1):
                squares[i][j].unload()
                for chn in channels[i][j]:
                    chn.unload()
                    chn.p_zoom = -1

    #SYMBOLS
    pause_symbol = 32  # SPACE
    jump_forward_symbol = 116  # T
    toggle_display_symbol = 121  # Y
    toggle_subdisplay_symbol = 117  # U
    clear_infowindows_symbol = 105  # I
    SETTINGS_KEY = 109  # M
    enter_symbol = 65293  # ENTER
    #endSYMBOLS

    #KEYS
    up_key = pyglet.window.key.UP
    down_key = pyglet.window.key.DOWN
    right_key = pyglet.window.key.RIGHT
    left_key = pyglet.window.key.LEFT
    zoom_key = pyglet.window.key.R
    dezoom_key = pyglet.window.key.E
    sts_plus_key = pyglet.window.key.D
    sts_minus_key = pyglet.window.key.A
    fps_minus_key = pyglet.window.key.Z
    fps_plus_key = pyglet.window.key.X
    MODIFIER_CTRL = pyglet.window.key.P
    MODIFIER_SHIFT = pyglet.window.key.W
    MODIFIER_ALT = pyglet.window.key.MOD_ALT
    #endKEYS

    @window.event
    def on_key_press(symbol, modifiers):
        global paused, CURRENT_FRAME, DISPLAY_MODE, SETTINGS_MODE
        if symbol == SETTINGS_KEY:
            SETTINGS_MODE = not SETTINGS_MODE
        if SETTINGS_MODE:
            return
        focus = None

        if timebar.SELECTING and (48 <= symbol <= 57 or symbol == 65288):
            if symbol == 65288:
                if len(timebar.sel_string) > 0:
                    timebar.sel_string = timebar.sel_string[0:-1]
            else:
                n_sym = symbol - 48
                timebar.sel_string += str(n_sym)

        if focus is None:
            if symbol == pause_symbol:
                paused = not paused
            elif symbol == enter_symbol:
                if timebar.SELECTING:
                    timebar.SELECTING = False
                    CURRENT_FRAME = max(0, min(N_FRAMES, int(timebar.sel_string) - STEPS_TO_SEE_MAX))
                    channel_wipe()
                    log_search()
                    r_update(STEPS_TO_SEE_MAX)
                    channel_relight()
                    timebar.sel_string = ""
            elif symbol == toggle_display_symbol:
                DISPLAY_MODE += 1
                DISPLAY_MODE %= DISPLAY_MODE_MOD
                full_unload()
            elif symbol == toggle_subdisplay_symbol:
                DISPLAY_MODE_SND[DISPLAY_MODE] += 1
                DISPLAY_MODE_SND[DISPLAY_MODE] %= DISPLAY_MODE_SND_MOD[DISPLAY_MODE]
            elif symbol == clear_infowindows_symbol:
                info.clear()
                info.LOCKED = False
                info_list.clear()
                gc.collect()

        print(symbol)

    btn_down = {1: False, 2: False, 3: False, 4: False, 5: False}

    @window.event
    def on_mouse_release(x, y, button, modifiers):
        if SETTINGS_MODE:
            return
        print('m_release', button, x, y)
        btn_down[button] = False

    @window.event
    def on_mouse_press(x, y, button, modifiers):
        global CURRENT_FRAME
        if SETTINGS_MODE:
            return
        print('m_press', button, x, y)
        btn_down[button] = True
        focus = None

        has_secondary_info = False
        secondary_window = None
        for inf in info_list:
            if inf.y - inf.h < y < inf.y and inf.x - inf.w < x < inf.x:
                has_secondary_info = True
                secondary_window = inf

        if (window.p_width//2) - (timebar.w//2) < x < (window.p_width//2) + (timebar.w//2) and timebar.y - 6 < y < timebar.y + 6:
            focus = 0

            t_x = x - ((window.p_width//2) - (timebar.w//2))
            frame_scale = t_x/timebar.w
            new_frame = int(frame_scale * N_FRAMES)

            CURRENT_FRAME = max(0, min(N_FRAMES, new_frame - STEPS_TO_SEE_MAX))
            channel_wipe()
            log_search()
            r_update(STEPS_TO_SEE_MAX)
            channel_relight()

        elif timebar.label.y < y < timebar.label.y + timebar.label.content_height and timebar.label.x - timebar.label.content_width/2 < x < timebar.label.x + timebar.label.content_width/2:
            focus = 1
            timebar.sel_string = ""
            timebar.SELECTING = True

        elif info.y - info.h < y < info.y and info.x - info.w < x < info.x:
            focus = 2

        elif has_secondary_info:
            focus = 3

        info.unselect()
        for inf in info_list:
            inf.unselect()

        if focus == 2:
            info.select(x, y)

        if focus == 3:
            secondary_window.select(x, y)

        if (not (focus == 2 or focus == 3)) or (focus == 2 and not info.LOCKED):
            if keys[MODIFIER_CTRL] and (not (focus == 2 or focus == 3)):
                n_win = InfoWindow(window.p_width - 20, window.p_height - 20, 200, 300, "Info", 40)
                info_list.append(n_win)
                find_nearest_object(x, y, 100, False, snd_win=n_win)
                n_win.LOCKED = True
            else:
                find_nearest_object(x, y, 100, False)
                if focus is None or (focus == 2 and not info.LOCKED):
                    if not (info.channel is None and info.event is None):
                        info.LOCKED = True
                    else:
                        info.LOCKED = False

        print(focus, 'fcs')
        print(button)

    @window.event
    def on_resize(width, height):
        global half_w, half_h

        p_p_width = window.p_width
        p_p_height = window.p_height

        window.p_width = width
        window.p_height = height

        dx = window.p_width - p_p_width
        dy = window.p_height - p_p_height

        half_w = window.p_width//2
        half_h = window.p_height//2

        info.adjust(dx, dy)
        for inf in info_list:
            inf.adjust(dx, dy)

        pv_x = int((-ox - (half_w/zoom)) // B_WIDTH) - 1
        pv_x_2 = int(pv_x + ((2 * half_w/zoom) // B_WIDTH)) + 1
        pv_x = max(0, min(WIDTH + 1, pv_x))
        pv_x_2 = max(0, min(WIDTH + 1, pv_x_2))

        pv_y = int((-oy - (half_h/zoom)) // B_HEIGHT)
        pv_y_2 = int(pv_y + ((2 * half_h/zoom) // B_HEIGHT)) + 2
        pv_y = max(0, min(HEIGHT + 1, pv_y))
        pv_y_2 = max(0, min(HEIGHT + 1, pv_y_2))

        for i in range(pv_x, pv_x_2 + 1, 1):
            for j in range(pv_y, pv_y_2 + 1, 1):
                for chn in channels[i][j]:
                    chn.up_zoom()

        if SETTINGS_MODE:
            settings_menu_adjust()
            return

    u_border = (HEIGHT + 2) * B_HEIGHT
    d_border = 0
    r_border = int((WIDTH + 2) * B_WIDTH)
    l_border = 0

    def fire_event(event):
        for e in event.data:
            channels_linear[e[0]].trigger_spike(e[1])
            channels_linear[e[0]].add_event(e[1], event)
        try:
            squares[int((event.get_x())/B_WIDTH)][int((event.get_y())/B_HEIGHT)].trigger()
        except IndexError:
            print("Index error: event triggered square out of range.")
        except OverflowError:
            print("Overflow error: Squares access overflowed to numpy.infty.")

    event_rad = 15
    spike_rad = 10

    def draw_event(event):
        pyglet.gl.glColor3f(BG_R, BG_G, BG_B)
        oex = ((event.get_x()) + ox) * zoom + half_w
        oey = ((event.get_y()) + oy) * zoom + half_h
        d_main = False
        for d in event.data:
            if DISPLAY_MODE == 0:
                if STEPS_TO_SEE_MAX - STEPS_TO_SEE <= STEPS_TO_SEE_MAX - CURRENT_FRAME + d[1] < STEPS_TO_SEE_MAX:

                    if not d_main:
                        if DISPLAY_MODE_SND[DISPLAY_MODE] == 1:
                            index_ = STEPS_TO_SEE - CURRENT_FRAME + d[1]
                            if 0 <= index_ <= STEPS_TO_SEE//4:
                                brt = index_/(STEPS_TO_SEE//4)
                                pyglet.gl.glColor4f(BG_R, BG_G, BG_B, brt)
                        circle_fill(event.get_x(), event.get_y(), event_rad + 2)
                        pyglet.gl.glColor3f(event.r, event.g, event.b)
                        if DISPLAY_MODE_SND[DISPLAY_MODE] == 1:
                            index_ = STEPS_TO_SEE - CURRENT_FRAME + d[1]
                            if 0 <= index_ <= STEPS_TO_SEE:
                                n_index_ = int((index_/STEPS_TO_SEE)*(len(brightness_bank_2) - 1))
                                _ = brightness_bank_2[n_index_]
                                pyglet.gl.glColor4f(_[0], _[1], _[2], _[3])
                            else:
                                pyglet.gl.glColor4f(0, 0, 0, 0)
                        circle_fill(event.get_x(), event.get_y(), event_rad)

                    d_main = True

                    chn = channels_linear[d[0]]

                    prv_x = half_w + int((chn.x + ox) * zoom)
                    prv_y = half_h + int((chn.y + oy) * zoom)

                    if prv_x < 0 or prv_x > window.p_width or prv_y < 0 or prv_y > window.p_height:
                        if chn.p_ox != ox:
                            chn.adjust()

                    x = channels_linear[d[0]].x_queue[STEPS_TO_SEE_MAX - CURRENT_FRAME + d[1]]
                    y = channels_linear[d[0]].queue[STEPS_TO_SEE_MAX - CURRENT_FRAME + d[1]]

                    if prv_x < 0 or prv_x > window.p_width or prv_y < 0 or prv_y > window.p_height:
                        #x = prv_x
                        y = prv_y

                    circle_og(x, y, spike_rad)

                    vec_x = x - oex
                    vec_y = y - oey
                    vec_len = math.sqrt((vec_x * vec_x) + (vec_y * vec_y))
                    vec_x /= vec_len
                    vec_y /= vec_len

                    e_x = vec_x * event_rad * zoom
                    spike_x = vec_x * -spike_rad * zoom
                    e_y = vec_y * event_rad * zoom
                    spike_y = vec_y * -spike_rad * zoom

                    pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                    pyglet.gl.glVertex2f(oex + e_x, oey + e_y)
                    pyglet.gl.glVertex2f(x + spike_x, y + spike_y)
                    pyglet.gl.glEnd()

            if DISPLAY_MODE == 3:
                if STEPS_TO_SEE_MAX - CURRENT_FRAME + d[1] < STEPS_TO_SEE_MAX:

                    if not d_main:
                        circle_fill(event.get_x(), event.get_y(), event_rad + 2)
                        pyglet.gl.glColor3f(event.r, event.g, event.b)
                        if DISPLAY_MODE_SND[DISPLAY_MODE] == 1:
                            index_ = STEPS_TO_SEE_MAX - CURRENT_FRAME + d[1]
                            if 0 <= index_ <= STEPS_TO_SEE_MAX:
                                n_index_ = int((index_/STEPS_TO_SEE_MAX)*(len(brightness_bank_2) - 1))
                                _ = brightness_bank_2[n_index_]
                                pyglet.gl.glColor4f(_[0], _[1], _[2], _[3])
                            else:
                                pyglet.gl.glColor4f(0, 0, 0, 0)
                        circle_fill(event.get_x(), event.get_y(), event_rad)

                    d_main = True

                    chn = channels_linear[d[0]]

                    prv_x = half_w + int((chn.x + ox) * zoom)
                    prv_y = half_h + int((chn.y + oy) * zoom)

                    if prv_x < 0 or prv_x > window.p_width or prv_y < 0 or prv_y > window.p_height:
                        if chn.p_ox != ox:
                            chn.adjust()

                    x = prv_x
                    y = prv_y

                    circle_og(x, y, spike_rad)

                    vec_x = x - oex
                    vec_y = y - oey
                    vec_len = math.sqrt((vec_x * vec_x) + (vec_y * vec_y))
                    vec_x /= vec_len
                    vec_y /= vec_len

                    e_x = vec_x * event_rad * zoom
                    spike_x = vec_x * -spike_rad * zoom
                    e_y = vec_y * event_rad * zoom
                    spike_y = vec_y * -spike_rad * zoom

                    pyglet.gl.glBegin(pyglet.gl.GL_LINES)
                    pyglet.gl.glVertex2f(oex + e_x, oey + e_y)
                    pyglet.gl.glVertex2f(x + spike_x, y + spike_y)
                    pyglet.gl.glEnd()

    def handle_HS2():
        global event_index
        try:
            event = process_post_event(event_index, mmap_indices, mmap_xy, mmap_data)
        except Exception as e:
            print(str(e))
            return
        if THEME == 0:
            bank = 1.5
            event.r = 0.2 + random.uniform(0, 0.8)
            bank -= event.r - 0.2
            event.g = 0.2 + random.uniform(0, min(0.8, bank))
            bank -= event.g - 0.2
            event.b = 0.2 + bank
        else:
            bank = 1
            event.r = 0.14 + random.uniform(0, 0.4)
            bank -= event.r - 0.14
            event.g = 0.14 + random.uniform(0, min(0.4, bank))
            bank -= event.g - 0.14
            event.b = 0.14 + bank

        while event.min_frame_v <= CURRENT_FRAME + 10:
            fire_event(event)
            event_index += 1
            try:
                event = process_post_event(event_index, mmap_indices, mmap_xy, mmap_data)
            except Exception as e:
                print(str(e))
                return
            if THEME == 0:
                bank = 1.5
                event.r = 0.2 + random.uniform(0, 0.8)
                bank -= event.r - 0.2
                event.g = 0.2 + random.uniform(0, min(0.8, bank))
                bank -= event.g - 0.2
                event.b = 0.2 + bank
            else:
                bank = 0.75
                event.r = 0.14 + random.uniform(0, 0.4)
                bank -= event.r - 0.14
                event.g = 0.14 + random.uniform(0, min(0.4, bank))
                bank -= event.g - 0.14
                event.b = 0.14 + bank

    keys = pyglet.window.key.KeyStateHandler()
    window.push_handlers(keys)

    prev_time = time.time()
    dcount = 0
    fps = "---"

    target_fps = 120
    actual_fps = 1

    def drawing_routine(x1b, y1b, x2b, y2b):
        for i in range(x1b, x2b + 1, 1):
            for j in range(y1b, y2b + 1, 1):
                for k in channels[i][j]:
                    for l in k.events:
                        if l[0] <= CURRENT_FRAME:
                            event = l[1]
                            if not event.TAGGED:
                                if DISPLAY_MODE == 0 or DISPLAY_MODE == 3:
                                    draw_event(event)
                                event.TAGGED = True

    def loading_routine(x1, y1, x2, y2, x1b, y1b, x2b, y2b, y_shift, zoom_shift):

        for i in range(x1, x2 + 1, 1):
            for j in range(y1, y2 + 1, 1):
                if not (x1b <= i <= x2b and y1b <= j <= y2b):
                    for k in channels[i][j]:
                        k.unload()

        for i in range(x1b, x2b + 1, 1):
            for j in range(y1b, y2b + 1, 1):
                if not (x1 <= i <= x2 and y1 <= j <= y2):
                    for k in channels[i][j]:
                        k.load()
                elif zoom_shift and (x1 <= i <= x2 and y1 <= j <= y2):
                    for k in channels[i][j]:
                        k.zoom_adjust()
                elif y_shift and (x1 <= i <= x2 and y1 <= j <= y2):
                    for k in channels[i][j]:
                        k.y_adjust()

    DISPLAY_MODE = 0
    DISPLAY_MODE_SND = []
    DISPLAY_MODE_MOD = 4
    for i in range(DISPLAY_MODE_MOD):
        DISPLAY_MODE_SND.append(0)
    DISPLAY_MODE_SND_MOD = [2, 1, len(mmap_totals_proc), 2]
    label_fontsize = 12
    max_tFPS = 244

    label = pyglet.text.Label("Fps: " + fps + (" (paused)" if paused else ""),
                              font_name='Lucida Console',
                              font_size=label_fontsize,
                              x=10, y=10,
                              color=(int(BASE_R * 255), int(BASE_G * 255), int(BASE_B * 255), 255),
                              anchor_x='left', anchor_y='bottom')

    label2 = pyglet.text.Label("TFps: " + str(target_fps) + "/" + str(max_tFPS),
                               font_name='Lucida Console',
                               font_size=label_fontsize,
                               x=10, y=10 + (2 * label_fontsize),
                               color=(int(BASE_R * 255), int(BASE_G * 255), int(BASE_B * 255), 255),
                               anchor_x='left', anchor_y='bottom')

    label3 = pyglet.text.Label("Zoom: " + str(round(zoom, 3)),
                               font_name='Lucida Console',
                               font_size=label_fontsize,
                               x=10, y=10 + (4 * label_fontsize),
                               color=(int(BASE_R * 255), int(BASE_G * 255), int(BASE_B * 255), 255),
                               anchor_x='left', anchor_y='bottom')

    label4 = pyglet.text.Label("STS: " + str(STEPS_TO_SEE) + '/' + str(STEPS_TO_SEE_MAX),
                               font_name='Lucida Console',
                               font_size=label_fontsize,
                               x=10, y=10 + (6 * label_fontsize),
                               color=(int(BASE_R * 255), int(BASE_G * 255), int(BASE_B * 255), 255),
                               anchor_x='left', anchor_y='bottom')

    label5 = pyglet.text.Label("Display: " + str(DISPLAY_MODE + 1) + "/" + str(DISPLAY_MODE_MOD) + " (" + str(DISPLAY_MODE_SND[DISPLAY_MODE] + 1) + "/" + str(DISPLAY_MODE_SND_MOD[DISPLAY_MODE]) + ")",
                               font_name='Lucida Console',
                               font_size=label_fontsize,
                               x=10, y=10 + (8 * label_fontsize),
                               color=(int(BASE_R * 255), int(BASE_G * 255), int(BASE_B * 255), 255),
                               anchor_x='left', anchor_y='bottom')

    pv_x = int((-ox - (half_w / zoom)) // B_WIDTH) - 1
    pv_x_2 = int(pv_x + ((2 * half_w / zoom) // B_WIDTH)) + 1
    pv_x = max(0, min(WIDTH + 1, pv_x))
    pv_x_2 = max(0, min(WIDTH + 1, pv_x_2))

    pv_y = int((-oy - (half_h / zoom)) // B_HEIGHT)
    pv_y_2 = int(pv_y + ((2 * half_h / zoom) // B_HEIGHT)) + 2
    pv_y = max(0, min(HEIGHT + 1, pv_y))
    pv_y_2 = max(0, min(HEIGHT + 1, pv_y_2))

    for i in range(pv_x, pv_x_2 + 1, 1):
        for j in range(pv_y, pv_y_2 + 1, 1):
            for chn in channels[i][j]:
                chn.load()

    def widget_jump(x, y):
        global ox, oy
        ox = -x
        oy = -y

    prev_ox = ox
    prev_oy = oy

    def update(dt):
        global prev_ox, prev_oy, prev_time, dcount, fps, mod_count, AF_DATA, dx, dy, ox, oy, event_index, dzoom, CURRENT_FRAME, half_w, half_h, label_fontsize, STEPS_TO_SEE, target_fps, MAIN_MEMORY_DATA, MAIN_MEMORY_COUNT, MAIN_MEMORY_COUNT_PLUS
        dtime = time.time() - prev_time
        dcount += 1
        if dtime >= 0.5:
            fps = str(round(dcount/dtime, 2))
            dcount = 0
            prev_time = time.time()
        window.dispatch_events()
        pyglet.gl.glClearColor(BG_R, BG_G, BG_B, 1)
        window.clear()

        if not paused:
            if MEMORY_MODEL == 0:
                data = probe.Read(CURRENT_FRAME, CURRENT_FRAME + 1)
            elif MEMORY_MODEL == 1:
                data = MAIN_MEMORY_DATA[MAIN_MEMORY_COUNT]
                MAIN_MEMORY_COUNT += 1
                if MAIN_MEMORY_COUNT == MAIN_MEMORY_WINDOW:
                    MAIN_MEMORY_COUNT = 0
                    MAIN_MEMORY_COUNT_PLUS += 1
                    MAIN_MEMORY_DATA = probe.Read(MAIN_MEMORY_COUNT_PLUS * MAIN_MEMORY_WINDOW, (1 + MAIN_MEMORY_COUNT_PLUS) * MAIN_MEMORY_WINDOW)
            AF_DATA[mod_count] = data
            mod_count += 1
            mod_count %= STEPS_TO_SEE_MAX

        if keys[dezoom_key]:
            if keys[MODIFIER_SHIFT]:
                label_fontsize -= 1
                label_fontsize = max(1, label_fontsize)
            else:
                dzoom -= 0.001 + abs(dzoom) * 0.1
                dzoom = max(dzoom, -0.05)
        elif keys[zoom_key]:
            if keys[MODIFIER_SHIFT]:
                label_fontsize += 1
                label_fontsize = max(1, label_fontsize)
            else:
                dzoom += 0.001 + abs(dzoom) * 0.1
                dzoom = max(dzoom, 0.05)
        else:
            dzoom -= (dzoom * 0.2)
            if abs(dzoom) < 0.04:
                dzoom = 0

        if keys[up_key] or keys[down_key]:
            if keys[up_key]:
                dy += abs(dy * 0.1) + 1
                dy = min(dy, 100)
            elif keys[down_key]:
                dy -= abs(dy * 0.1) + 1
                dy = max(dy, -100)
        else:
            dy -= (dy * 0.2)
            if abs(dy) < 1:
                dy = 0

        if keys[right_key] or keys[left_key]:
            if keys[left_key]:
                dx += abs(dx * 0.1) + 1
                dx = min(dx, 100)
            elif keys[right_key]:
                dx -= abs(dx * 0.1) + 1
                dx = max(dx, -100)
        else:
            dx -= (dx * 0.2)
            if abs(dx) < 1:
                dx = 0

        if keys[sts_minus_key] or keys[sts_plus_key]:
            if keys[sts_minus_key]:
                STEPS_TO_SEE -= 1
                STEPS_TO_SEE = max(1, STEPS_TO_SEE)

            if keys[sts_plus_key]:
                STEPS_TO_SEE += 1
                STEPS_TO_SEE = min(STEPS_TO_SEE_MAX, STEPS_TO_SEE)

            pv_x = int((-prev_ox - (half_w / zoom)) // B_WIDTH) - 1
            pv_x_2 = int(pv_x + ((2 * half_w / zoom) // B_WIDTH)) + 1
            pv_x = max(0, min(WIDTH + 1, pv_x))
            pv_x_2 = max(0, min(WIDTH + 1, pv_x_2))

            pv_y = int((-prev_oy - (half_h / zoom)) // B_HEIGHT)
            pv_y_2 = int(pv_y + ((2 * half_h / zoom) // B_HEIGHT)) + 2
            pv_y = max(0, min(HEIGHT + 1, pv_y))
            pv_y_2 = max(0, min(HEIGHT + 1, pv_y_2))

            for i in range(pv_x, pv_x_2 + 1, 1):
                for j in range(pv_y, pv_y_2 + 1, 1):
                    for chn in channels[i][j]:
                        chn.adjust()

        if keys[fps_minus_key] or keys[fps_plus_key]:
            if keys[fps_minus_key]:
                target_fps -= 1
                target_fps = max(0, target_fps)

            if keys[fps_plus_key]:
                target_fps += 1
                target_fps = min(max_tFPS, target_fps)

        i_pv_x = int((-prev_ox - (half_w/zoom)) // B_WIDTH) - 1
        i_pv_x_2 = int(i_pv_x + ((2 * half_w/zoom) // B_WIDTH)) + 1
        i_pv_x = max(0, min(WIDTH + 1, i_pv_x))
        i_pv_x_2 = max(0, min(WIDTH + 1, i_pv_x_2))

        i_pv_y = int((-prev_oy - (half_h/zoom)) // B_HEIGHT)
        i_pv_y_2 = int(i_pv_y + ((2 * half_h/zoom) // B_HEIGHT)) + 2
        i_pv_y = max(0, min(HEIGHT + 1, i_pv_y))
        i_pv_y_2 = max(0, min(HEIGHT + 1, i_pv_y_2))

        x_shift = False
        y_shift = False
        zoom_shift = False

        if dx != 0:
            window.scroll_x(int(dx))
            x_shift = True

        if dy != 0:
            window.scroll_y(int(dy))
            y_shift = True

        if dzoom != 0:
            handle_zoom(1 + dzoom)
            zoom_shift = True

        pv_x = int((-ox - (half_w/zoom)) // B_WIDTH) - 1
        pv_x_2 = int(pv_x + ((2 * half_w/zoom) // B_WIDTH)) + 1
        pv_x = max(0, min(WIDTH + 1, pv_x))
        pv_x_2 = max(0, min(WIDTH + 1, pv_x_2))

        pv_y = int((-oy - (half_h/zoom)) // B_HEIGHT)
        pv_y_2 = int(pv_y + ((2 * half_h/zoom) // B_HEIGHT)) + 2
        pv_y = max(0, min(HEIGHT + 1, pv_y))
        pv_y_2 = max(0, min(HEIGHT + 1, pv_y_2))

        loading_routine(i_pv_x, i_pv_y, i_pv_x_2, i_pv_y_2, pv_x, pv_y, pv_x_2, pv_y_2, y_shift, zoom_shift)

        # bbox it

        if DISPLAY_MODE == 1:
            for i in range(pv_x, pv_x_2 + 1, 1):
                for j in range(pv_y, pv_y_2 + 1, 1):
                    squares[i][j].update()
            for i in range(pv_x, pv_x_2 + 1, 1):
                for j in range(pv_y, pv_y_2 + 1, 1):
                    squares[i][j].draw()

        if DISPLAY_MODE == 2:
            for i in range(pv_x, pv_x_2 + 1, 1):
                for j in range(pv_y, pv_y_2 + 1, 1):
                    squares[i][j].brightness = mmap_totals_proc[DISPLAY_MODE_SND[DISPLAY_MODE]][i][j]
            for i in range(pv_x, pv_x_2 + 1, 1):
                for j in range(pv_y, pv_y_2 + 1, 1):
                    squares[i][j].draw()

        pyglet.gl.glColor3f(BASE_R, BASE_G, BASE_B)

        if -window.p_width//2 < int((l_border + ox) * zoom) < window.p_width//2:
            pyglet.gl.glBegin(pyglet.gl.GL_LINES)
            pyglet.gl.glVertex2f(half_w + int((l_border + ox) * zoom), 0)
            pyglet.gl.glVertex2f(half_w + int((l_border + ox) * zoom), window.p_height)
            pyglet.gl.glEnd()

        if -window.p_width//2 < int((r_border + ox) * zoom) < window.p_width//2:
            pyglet.gl.glBegin(pyglet.gl.GL_LINES)
            pyglet.gl.glVertex2f(half_w + int((r_border + ox) * zoom), 0)
            pyglet.gl.glVertex2f(half_w + int((r_border + ox) * zoom), window.p_height)
            pyglet.gl.glEnd()

        if -window.p_height//2 < int((d_border + oy) * zoom) < window.p_height//2:
            pyglet.gl.glBegin(pyglet.gl.GL_LINES)
            pyglet.gl.glVertex2f(0, half_h + int((d_border + oy) * zoom))
            pyglet.gl.glVertex2f(window.p_width, half_h + int((d_border + oy) * zoom))
            pyglet.gl.glEnd()

        if -window.p_height//2 < int((u_border + oy) * zoom) < window.p_height//2:
            pyglet.gl.glBegin(pyglet.gl.GL_LINES)
            pyglet.gl.glVertex2f(0, half_h + int((u_border + oy) * zoom))
            pyglet.gl.glVertex2f(window.p_width, half_h + int((u_border + oy) * zoom))
            pyglet.gl.glEnd()

        if DISPLAY_MODE == 0 or DISPLAY_MODE == 3:

            for i in range(pv_x, pv_x_2 + 1, 1):
                for j in range(pv_y, pv_y_2 + 1, 1):
                    for chn in channels[i][j]:
                        if not paused:
                            chn.update(data[chn.index])
                        else:
                            chn.recalc()

            drawing_routine(pv_x, pv_y, pv_x_2, pv_y_2)

            if DISPLAY_MODE == 0:

                if DISPLAY_MODE_SND[DISPLAY_MODE] == 1:
                    for i in range(pv_x, pv_x_2 + 1, 1):
                        for j in range(pv_y, pv_y_2 + 1, 1):
                            for chn in channels[i][j]:
                                if chn.l_total == 0:
                                    window.draw_multiline_decay(chn.x_queue, chn.queue)
                                else:
                                    window.draw_multiline_2_decay(chn.x_queue, chn.queue, chn.l_queue)
                else:
                    for i in range(pv_x, pv_x_2 + 1, 1):
                        for j in range(pv_y, pv_y_2 + 1, 1):
                            for chn in channels[i][j]:
                                if chn.l_total == 0:
                                    window.draw_multiline(chn.x_queue, chn.queue)
                                else:
                                    window.draw_multiline_2(chn.x_queue, chn.queue, chn.l_queue)
                if render_mode == 0:
                    window.batch_draw()

        if not paused:
            handle_HS2()
            CURRENT_FRAME = CURRENT_FRAME + 1
            CURRENT_FRAME = min(N_FRAMES - 1, max(0, CURRENT_FRAME))

        label.x = 10
        label.y = 10
        label.font_size = label_fontsize
        label.text = "Fps: " + fps + (" (paused)" if paused else "")
        label.draw()

        label2.x = 10
        label2.y = 10 + (2 * label_fontsize)
        label2.font_size = label_fontsize
        label2.text = "TFps: " + str(target_fps) + "/" + str(max_tFPS)
        label2.draw()

        label3.x = 10
        label3.y = 10 + (4 * label_fontsize)
        label3.font_size = label_fontsize
        label3.text = "Zoom: " + str(round(zoom, 3))
        label3.draw()

        label4.x = 10
        label4.y = 10 + (6 * label_fontsize)
        label4.font_size = label_fontsize
        label4.text = "STS: " + str(STEPS_TO_SEE) + '/' + str(STEPS_TO_SEE_MAX)
        label4.draw()

        label5.x = 10
        label5.y = 10 + (8 * label_fontsize)
        label5.font_size = label_fontsize
        label5.text = "Display: " + str(DISPLAY_MODE + 1) + "/" + str(DISPLAY_MODE_MOD) + " (" + str(DISPLAY_MODE_SND[DISPLAY_MODE] + 1) + "/" + str(DISPLAY_MODE_SND_MOD[DISPLAY_MODE]) + ")"
        label5.draw()

        info.draw()
        for inf in info_list:
            inf.draw()

        timebar.update(CURRENT_FRAME - 1, max(label.content_width, label2.content_width, label3.content_width))
        timebar.draw()

        pyglet.gl.glFlush()
        # gc.collect()  # Cleans memory leaks (heavy performance cost)
        prev_ox = ox
        prev_oy = oy
        if not SETTINGS_MODE:
            pyglet.clock.schedule_once(update, 1 / target_fps)
        else:
            pyglet.clock.schedule_once(settings_update, 1 / 60)

    #This update method only called for re-adjustement
    def r_update(numtimes):
        global prev_time, dcount, fps, mod_count, AF_DATA, dx, dy, ox, oy, event_index, dzoom, CURRENT_FRAME, half_w, half_h, label_fontsize, STEPS_TO_SEE, target_fps

        for _ in range(numtimes + 1):
            data = probe.Read(CURRENT_FRAME, CURRENT_FRAME + 1)
            AF_DATA[mod_count] = data
            mod_count += 1
            mod_count %= STEPS_TO_SEE_MAX

            handle_HS2()
            CURRENT_FRAME = CURRENT_FRAME + 1

        # bbox it

        pv_x = int((-ox - (half_w / zoom)) // B_WIDTH) - 1
        pv_x_2 = int(pv_x + ((2 * half_w / zoom) // B_WIDTH)) + 1
        pv_x = max(0, min(WIDTH + 1, pv_x))
        pv_x_2 = max(0, min(WIDTH + 1, pv_x_2))

        pv_y = int((-oy - (half_h / zoom)) // B_HEIGHT)
        pv_y_2 = int(pv_y + ((2 * half_h / zoom) // B_HEIGHT)) + 2
        pv_y = max(0, min(HEIGHT + 1, pv_y))
        pv_y_2 = max(0, min(HEIGHT + 1, pv_y_2))

        for i in range(pv_x, pv_x_2 + 1, 1):
            for j in range(pv_y, pv_y_2 + 1, 1):
                for chn in channels[i][j]:
                    chn.p_zoom = -1

        dcount = 0
        prev_time = time.time()

        # gc.collect()  # Cleans memory leaks (heavy performance cost)

    class KeySet:

        def __init__(self, k, text, upper):
            self.k = k
            self.l1 = self.label = pyglet.text.Label(
                        text + " ",
                        font_name='Lucida Console',
                        font_size=14,
                        x=window.p_width//2, y=upper + d_set_label,
                        color=(int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), 255),
                        anchor_x='right',
                        anchor_y='bottom')
            self.l2 = self.label = pyglet.text.Label(
                        " " + pyglet.window.key.symbol_string(k),
                        font_name='Lucida Console',
                        font_size=14,
                        x=window.p_width//2, y=upper + d_set_label,
                        color=(int(BASE_R*255), int(BASE_G*255), int(BASE_B*255), 255),
                        anchor_x='left',
                        anchor_y='bottom')

        def update(self):
            self.l1.x = window.p_width // 2
            self.l2.x = window.p_width // 2

        def draw(self):
            self.l1.draw()
            self.l2.draw()

    key_labels = []
    _mset = 0
    d_set_label = 3

    def key_draw(k, s):
        global key_labels, _mset
        ks = KeySet(k, s, _mset)
        key_labels.append(ks)
        _mset = max(ks.l1.y + ks.l1.content_height, ks.l2.y + ks.l2.content_height) + d_set_label

    key_draw(pause_symbol, "Pause")
    # key_draw(jump_forward_symbol, "")
    key_draw(toggle_display_symbol, "Toggle Display")
    key_draw(toggle_subdisplay_symbol, "Toggle Subdisplay")
    key_draw(clear_infowindows_symbol, "Clear Info Windows")
    key_draw(enter_symbol, "Confirm")

    key_draw(up_key, "Scroll Up")
    key_draw(down_key, "Scroll Down")
    key_draw(right_key, "Scroll Right")
    key_draw(left_key, "Scroll Left")
    key_draw(zoom_key, "Zoom In")
    key_draw(dezoom_key, "Zoom Out")
    key_draw(sts_plus_key, "Increase Visible Steps")
    key_draw(sts_minus_key, "Decrease Visible Steps")
    key_draw(fps_minus_key, "Decrease Target FPS")
    key_draw(fps_plus_key, "Increase Target FPS")
    key_draw(MODIFIER_CTRL, "Modifier 1")
    key_draw(SETTINGS_KEY, "Settings Menu")
    # key_draw(MODIFIER_SHIFT, "")
    # key_draw(MODIFIER_ALT, "")

    setting_scroll = 0

    def settings_menu_adjust_2(fnt):
        for kl in key_labels:
            kl.l1.font_size = fnt
            kl.l2.font_size = fnt
        c_scroll = d_set_label
        for kl in key_labels:
            kl.l1.y = c_scroll
            kl.l2.y = c_scroll
            c_scroll += max(kl.l1.content_height, kl.l2.content_height) + d_set_label

    def settings_menu_adjust():
        fnt = key_labels[-1].l1.font_size
        while key_labels[-1].l1.y + key_labels[-1].l1.content_height > window.p_height:
            fnt -= 1
            fnt = max(fnt, 1)
            settings_menu_adjust_2(fnt)
            if fnt == 1:
                break
        while key_labels[-1].l1.y + key_labels[-1].l1.content_height < window.p_height:
            fnt += 1
            fnt = max(fnt, 1)
            settings_menu_adjust_2(fnt)
        fnt -= 1
        fnt = max(fnt, 1)
        settings_menu_adjust_2(fnt)

    def settings_update(dt):
        window.dispatch_events()
        pyglet.gl.glClearColor(BG_R, BG_G, BG_B, 1)
        window.clear()

        pyglet.gl.glBegin(pyglet.gl.GL_LINES)
        pyglet.gl.glColor4f(BASE_R, BASE_G, BASE_B, 1)
        pyglet.gl.glVertex2f((window.p_width // 2), 0)
        pyglet.gl.glVertex2f((window.p_width // 2), window.p_height)
        for k in key_labels:
            x1 = 0
            x2 = window.p_width
            y = k.l1.y - d_set_label
            pyglet.gl.glVertex2f(x1, y)
            pyglet.gl.glVertex2f(x2, y)
        pyglet.gl.glEnd()

        for k in key_labels:
            k.update()
            k.draw()

        pyglet.gl.glFlush()
        # gc.collect()  # Cleans memory leaks (heavy performance cost)
        if not SETTINGS_MODE:
            gc.collect()
            pyglet.clock.schedule_once(update, 1 / target_fps)
        else:
            pyglet.clock.schedule_once(settings_update, 1 / 60)

    pyglet.clock.schedule_once(update, 1/target_fps)

    pyglet.app.run()
