import chimera.extension

class VolumeCut_EMO(chimera.extension.EMO):
    def name(self):
        return 'VolumeCut'

    def description(self):
        return 'Cut Volume surrounding a protein chain.'

    def categories(self):
        return ['Utilities']

    def icon(self):
        return self.path('cutvol.tiff')

    def activate(self):
	# Call the 'mainchain' function in the "__init__.py" module.
        self.module('gui').show_vol_dialog()

chimera.extension.manager.registerExtension(VolumeCut_EMO(__file__))
