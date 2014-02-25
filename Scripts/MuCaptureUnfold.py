#!/usr/bin/env python
# ==============================================================================
#  Description:
#       A simple code to try my first unfolding
#
#  Author: Anthony Hillairet
#
# ==============================================================================



from ROOT import gROOT, gRandom, gStyle, TH1, TH1D, cout, TFile, TCanvas, TLegend, TGraph
from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
# from ROOT import RooUnfoldSvd
# from ROOT import RooUnfoldTUnfold
from optparse import OptionParser

def Openfile(inputopt, defaultTitle):
	inputSplit = inputopt.split(':')
	thefile = TFile(inputSplit[0])
	if thefile.IsZombie():
		print "ERROR: Problem when opening",inputSplit[0]
		sys.exit(1)
	if len(inputSplit) > 1:
		thetitle = ':'.join(inputSplit[1:])
	else:
		thetitle = defaultTitle

	return (thefile, thetitle)


class Unfolder:


	def __init__(self,options):

		fail = False
		required = ["unfold", "training", "output"]
		for opt in required:
			if not hasattr(options, opt):
				print "--" + opt + " option must be specified"
				fail = True
		
		if fail:
			sys.exit(1)

		# Open input files
		(self.unfold, self.unfoldTitle) = Openfile(options.unfold, "input")
		(self.training, self.trainingTitle) = Openfile(options.training, "training")

		# Open output file
		if options.output.find('.pdf') == -1:
		    self.output = options.output+'.pdf'
		else:
			self.output = options.output

		self.Canv = TCanvas()
		self.Canv.Print(self.output+'[')

		if self.unfold.Get("MuCapture/anDnLateResponse") == None:
			self.unfoldIsMC = False
		else:
			self.unfoldIsMC = True

		self.unfoldresult = None

		self.nbiterations = 10



	
	def PrintMeasured(self):

		def DrawMeasuredFromResponse(thefile, thetitle):
			response = thefile.Get("MuCapture/anDnLateResponse")
			hMeas = response.Hmeasured()
			if not thetitle == None:
				hMeas.SetTitle(thetitle+" measured spectrum")
			hMeas.Draw()
		
		DrawMeasuredFromResponse(self.training, self.trainingTitle)
		self.Canv.Print(self.output)
		self.Canv.Clear()

		if self.unfoldIsMC:
			DrawMeasuredFromResponse(self.unfold, self.unfoldTitle)
			self.Canv.Print(self.output)
			self.Canv.Clear()


	def PrintResponseMatrix(self):

		def DrawMeasuredFromResponse(thefile, thetitle):
			matrix = thefile.Get("MuCapture/LateResponse/MCMeasVsTruthMomentum")
			if not thetitle == None:
				matrix.SetTitle(thetitle+" response matrix")
			matrix.Draw("colz")
		
		DrawMeasuredFromResponse(self.training, self.trainingTitle)
		self.Canv.Print(self.output)
		self.Canv.Clear()



	def UnfoldBayes(self):
		response = self.training.Get("MuCapture/anDnLateResponse")

		# TODO: This is temporary until MeasuredMomentum is filled for both data and MC
		tmpresponse = self.unfold.Get("MuCapture/anDnLateResponse")
		unfoldMeas = tmpresponse.Hmeasured()

		self.unfoldresult = RooUnfoldBayes(response, unfoldMeas, self.nbiterations);
		self.unfoldname = "Bayes"


	def PrintResultsMC(self):
		if self.unfoldresult == None:
			print "ERROR: No unfolding performed yet !!!"
			sys.exit(1)

		hReco= self.unfoldresult.Hreco()

		if self.unfoldIsMC:
			unfoldTruth = self.unfold.Get("MuCapture/LateResponse/MCTruthMomentum")
			unfoldTruth.SetLineColor(2)
			leg = TLegend(0.65,0.5,0.98,0.7)

			MaxhReco = hReco.GetMaximum()
			unfoldTruth.Draw()
			hReco.Draw("same")
			leg.AddEntry(unfoldTruth, 'True input spectrum','L')
			leg.AddEntry(hReco, 'Unfolded spectrum', 'L')
			unfoldTruth.SetTitle("%s spectrum unfolded with %s spectrum using RooUnfold%s" % (self.unfoldTitle, self.trainingTitle, self.unfoldname))
			unfoldTruth.SetAxisRange(0.0,(1.1*MaxhReco),'Y')
			leg.Draw()
		else:
			hReco.Draw()
			hReco.SetTitle(self.trainingTitle+" spectrum unfolded using RooUnfold"+self.unfoldname)
		self.Canv.Update()
		self.Canv.Print(self.output)
		self.Canv.Clear()

	def MakeChi2VsIterationPlot(self):
		self.Chi2VsIter = TGraph()

	def AddChi2VsIterationPoint(self):
		if self.unfoldresult == None:
			print "ERROR: No unfolding performed yet !!!"
			sys.exit(1)

		if not self.unfoldIsMC:
			print "WARNING: There is no truth information to compare to. No Chi2 point added."
			return

		unfoldTruth = self.unfold.Get("MuCapture/LateResponse/MCTruthMomentum")
		chi2 = self.unfoldresult.Chi2(unfoldTruth)
		nbPt = self.Chi2VsIter.GetN()
		self.Chi2VsIter.SetPoint(nbPt, self.nbiterations, chi2)


	def DrawChi2VsIterationPlot(self):
		self.Canv.Clear()
		self.Chi2VsIter.SetTitle("chi^{2} between unfolded and true %s spectra, using %s spectrum;Number of iterations;chi^{2}" % (self.unfoldTitle,self.trainingTitle))
		self.Chi2VsIter.SetMarkerStyle(5)
		self.Chi2VsIter.SetMarkerColor(4)
		self.Chi2VsIter.SetLineColor(4)
		self.Chi2VsIter.Draw("APL")
		self.Canv.Update()
		self.Canv.Print(self.output)
		self.Canv.Clear()

	def Close(self):
		self.Canv.Print(self.output+"]")



# if __name__ == "__main__":


# OptionParser.format_description = lambda self, formatter: "\n".join(map(textwrap.fill, self.description.split("\n"))) + "\n"
# parser = OptionParser(description=description)
parser = OptionParser()
parser.add_option("-t", "--training", dest="training", help="Root file containing the rooUnfold response.")
parser.add_option("-u", "--unfold", dest="unfold", help="Root file containing the measured spectrum you want to unfold")
parser.add_option("-o", "--output", dest="output", help="Output file name for the pdf file. Don't put the suffix '.pdf'", default="output")
parser.add_option("-i", "--iterationcheck", dest="iterationcheck", help="Plot the unfolding results versus the number of iterations.", action="store_true")
(options, args) = parser.parse_args()

gROOT.SetBatch()
gStyle.SetOptStat(10)
gStyle.SetTitleW(1.0)

Unf = Unfolder(options)
# Study of unfolding results vs number of iterations
if options.iterationcheck:
	Unf.MakeChi2VsIterationPlot()
	for it in range(1,16):
		Unf.nbiterations = it
		Unf.UnfoldBayes()
		Unf.unfoldname = "Bayes, %d iterations" % it
		Unf.PrintResultsMC()
		Unf.AddChi2VsIterationPoint()
	Unf.DrawChi2VsIterationPlot()
# Default output
else:
	Unf.PrintMeasured()
	Unf.PrintResponseMatrix()
	Unf.UnfoldBayes()
	Unf.PrintResultsMC()
Unf.Close()
