from psychopy import core, visual, event, gui, monitors, data, misc


scnWidth, scnHeight = (1920, 1080)
screen_width        = 53.5 # centimeters 
screen_dist         = 55.0
white               = [255,255,255] # background screen color
grey                = [128,128,128]
black               = [0,0,0] # background screen color



mon = monitors.Monitor('behavlab', width=screen_width, distance=screen_dist)
mon.setSizePix((scnWidth, scnHeight))
win = visual.Window((scnWidth, scnHeight),color=white,colorSpace='rgb255',monitor=mon,fullscr= True ,units='deg')
win.update()


fixcross_txt = visual.TextStim(win, text='+', colorSpace = 'rgb255', color=[0,0,0], units = "deg", pos = [-10,0], height=0.8)
fixcross_png = visual.ImageStim(win, image = 'cross.png', units='deg', pos=[10,0], size=0.8)

cue_txt = visual.TextStim(win, text='*', colorSpace = 'rgb255', color=[0,0,0], units = "deg", pos = [-10,0], height=0.8)
cue_png = visual.ImageStim(win, image = 'asterisk_black.png', units='deg', pos=[10,0], size=0.8)


fixcross_txt.draw()
fixcross_png.draw()
win.flip()
event.waitKeys()


cue_txt.draw()
cue_png.draw()
win.flip()
event.waitKeys()



win.close()
