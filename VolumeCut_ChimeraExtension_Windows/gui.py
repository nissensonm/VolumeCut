import warnings
import chimera
from chimera.baseDialog import ModelessDialog
import VolumeViewer

# -----------------------------------------------------------------------------
#
class Vol_Cut_Dialog(ModelessDialog):

  title = 'Volume_Cut'
  name = 'Volume_Cut'
  buttons = ('Cut Volume', 'Close',)
  help = 'https://github.com/nissensonm/VolumeCut/'

    
# -----------------------------------------------------------------------------
#
  def fillInUI(self, parent):
    self.toplevel_widget = parent.winfo_toplevel()
    self.toplevel_widget.withdraw()
	
    parent.columnconfigure(0, weight = 1)

    from CGLtk import Hybrid
    import Pmw, Tkinter
	
    from chimera.widgets import MoleculeOptionMenu
    self.molMenu = MoleculeOptionMenu(parent, labelpos="w",
						label_text="Select PDB:")
    self.molMenu.grid(row = 1, column = 0, sticky = 'w')

    from chimera.widgets import MoleculeChainOptionMenu 
    self.chainMenu = MoleculeChainOptionMenu(parent, labelpos="w",
	                    label_text="Select Chain:")
    self.chainMenu.grid(row = 2, column = 0, sticky = 'w')
	
	
    from chimera.widgets import ModelOptionMenu 
    self.modelMenu = ModelOptionMenu(parent, labelpos="w",
	                    label_text="Select MRC:")
    self.modelMenu.grid(row = 3, column = 0, sticky = 'w')
	
    import Pmw
    self.radius = Pmw.EntryField(parent, labelpos="w",
                        label_text="Radius:")
    self.radius.grid(row = 4, column = 0, sticky = 'w')
	
    self.resolution = Pmw.EntryField(parent, labelpos="w",
                        label_text="Resolution:")
    self.resolution.grid(row = 5, column = 0, sticky = 'w')
	
# -----------------------------------------------------------------------------
#
  def CutVolume(self):
    #print ' mode: ' + str(self.mode.get()) 
    import VolCut as vc


    if self.molMenu.getvalue() is None or self.modelMenu.getvalue() is None \
                or self.chainMenu.getvalue() is None:
      warnings.warn('Must load valid PDB, MRC, or select a chain.')
      return
	  
    if type(self.modelMenu.getvalue()) is not VolumeViewer.volume.Volume:
      warnings.warn('Must select a valid MRC for the mrc')
      return
    
    try:
      float(self.radius.getvalue())
    except ValueError:
      warnings.warn('Radius must be a valid float value')
      return
	
    try:
      float(self.resolution.getvalue())
    except ValueError:
      warnings.warn('Resolution must be a valid float value')
      return
	
    vc.cutvol(self.molMenu.getvalue().id, self.modelMenu.getvalue().id, \
                float(self.radius.getvalue()), \
                float(self.resolution.getvalue()), \
                self.chainMenu.getvalue().chain)
   # if self.mode.get() == self.MEDIUM_RESOLUTION:
    #  print 'Cutting Low Resolution model'
    #  vc.cutvol(str(self.mm_pdb.variable.get()), str(self.mm_mrc.variable.get()), float(self.radius.getvalue()), float(self.resolution.getvalue()), "A")
  
# ---------------------------------------------------------------------------
#	  
	  
  # def Refresh(self):
    # print self.radius.getvalue() 
    # #print type(float(self.radius.getvalue()))
    # #print dir(self.molMenu)
    # #print self.molMenu.getvalue()
    # print dir(self.chainMenu.getvalue())
    # print self.chainMenu.getvalue()._name
    # print type(self.chainMenu.getvalue()._name)
    # print self.chainMenu.getvalue().chain
    # print self.chainMenu.getvalue().descriptiveName
    # #print dir(self.molMenu.getvalue())
    # #print self.molMenu.getvalue().id
	# #print self.chainMenu.getvalue()
    # #print self.chainMenu.getvar
    # #print self.chainMenu.info()
    # #print dir(self.chainMenu);
    # #import VolCut as vc
    # #active_PDB_models, active_MRC_models = vc.get_list_of_items()
	
    # #self.active_PDB_models = tuple([str(i) for i in active_PDB_models])
    # #self.active_MRC_models = tuple([str(i) for i in active_MRC_models])
    # #print(self.active_MRC_models);
    # #from CGLtk import Hybrid
    # #print 'eventually refresh the displayed PDB / MRC files.'
    # #self.mm_pdb.remove_all_entries(self)
    # #self.mm_pdb.add_entry(self, self.active_MRC_models)
    # #mm_pdb = Hybrid.Option_Menu(parent, 'PDB       ', *modes_PDB)
    # #mm_pdb.frame.grid(row = 1, column = 0, sticky = 'w')
    # #self.pdb_selected = mm_pdb.variable
    # #mm_pdb.add_callback(self.mode_changed_PDB)
	
  # ---------------------------------------------------------------------------
  #	  
def movement_mode_dialog(create = 0):
  from chimera import dialogs
  return dialogs.find(Vol_Cut_Dialog.name, create=create)	

  # ---------------------------------------------------------------------------
  #	  
def show_vol_dialog():
  from chimera import dialogs
  return dialogs.display(Vol_Cut_Dialog.name)
  
# -----------------------------------------------------------------------------

from chimera import dialogs
dialogs.register(Vol_Cut_Dialog.name,
                 Vol_Cut_Dialog, replace = 1)
