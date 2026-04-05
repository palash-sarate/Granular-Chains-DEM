import tkinter as tk
from typing import Callable, Optional

class PlaybackController:
    """Manages playback state, timer, and transport controls."""
    def __init__(self, root: tk.Misc, get_timesteps_cb: Callable, on_frame_change_cb: Callable):
        self.root = root
        self.get_timesteps_cb = get_timesteps_cb
        self.on_frame_change_cb = on_frame_change_cb
        
        self.playing = False
        self.job = None
        
        # UI References
        self.frame_slider = None
        self.fps_slider = None # This is now the FPS Entry
        self.play_btn = None
        self.frame_label = None
        self.loop_var = None

    def link_widgets(self, frame_slider, fps_slider, play_btn, frame_label, loop_var):
        self.frame_slider = frame_slider
        self.fps_slider = fps_slider
        self.play_btn = play_btn
        self.frame_label = frame_label
        self.loop_var = loop_var

    def update_status(self, idx: int):
        if not self.frame_label: return
        timesteps = self.get_timesteps_cb()
        total = len(timesteps)
        if total == 0:
            self.frame_label.config(text="No frames loaded")
            return
        ts = timesteps[idx]
        self.frame_label.config(text=f"Frame: {idx+1} / {total} (TS: {ts})")

    def jump_start(self):
        self.pause()
        if self.frame_slider and self.get_timesteps_cb():
            self.frame_slider.set(0)

    def jump_end(self):
        self.pause()
        ts = self.get_timesteps_cb()
        if self.frame_slider and ts:
            self.frame_slider.set(len(ts) - 1)

    def step_next(self):
        self.pause()
        if not self.frame_slider: return
        curr = int(self.frame_slider.get())
        ts = self.get_timesteps_cb()
        if curr < len(ts) - 1:
            self.frame_slider.set(curr + 1)
        elif self.loop_var and self.loop_var.get():
            self.frame_slider.set(0)

    def step_prev(self):
        self.pause()
        if not self.frame_slider: return
        curr = int(self.frame_slider.get())
        ts = self.get_timesteps_cb()
        if curr > 0:
            self.frame_slider.set(curr - 1)
        elif self.loop_var and self.loop_var.get():
            self.frame_slider.set(len(ts)-1)

    def pause(self):
        self.playing = False
        if self.job:
            self.root.after_cancel(self.job)
            self.job = None
        if self.play_btn:
            self.play_btn.config(text="▶ Play")

    def toggle_play(self):
        if self.playing:
            self.pause()
        else:
            if not self.get_timesteps_cb(): return
            self.playing = True
            if self.play_btn:
                self.play_btn.config(text="⏸ Pause")
            self._tick()

    def _tick(self):
        if not self.playing or not self.frame_slider: return
        
        curr = int(self.frame_slider.get())
        ts_list = self.get_timesteps_cb()
        total = len(ts_list)
        
        # Get FPS from text entry (default 10)
        try:
            fps_val = int(self.fps_slider.get())
        except (ValueError, AttributeError):
            fps_val = 10
            
        next_idx = curr + fps_val
        if next_idx >= total:
            if self.loop_var and self.loop_var.get():
                next_idx = 0
            else:
                self.pause()
                return

        self.frame_slider.set(next_idx)
        
        delay = 100 # 10 ticks per second
        self.job = self.root.after(delay, self._tick)
