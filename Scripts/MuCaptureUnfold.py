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

import sys

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
		required = ["unfold", "output"]
		for opt in required:
			if not hasattr(options, opt):
				print "--" + opt + " option must be specified"
				fail = True
		if options.training == "" and not options.multipletraining and not options.differencetraining:
			print "--training option must be specified"
			fail = True
		
		if fail:
			sys.exit(1)

		# Open input files
		(self.unfold, self.unfoldTitle) = Openfile(options.unfold, "input")
		if not options.training == "":
			self.InitTraining(options.training, "training")

		# Open output file
		if options.output.find('.pdf') == -1:
		    self.output = options.output+'.pdf'
		else:
			self.output = options.output

		if options.logy:
		    self.logy = True
		else:
		    self.logy = False

		self.Canv = TCanvas()
		self.Canv.Print(self.output+'[')

		self.ResponseName = options.responsename
		if self.unfold.Get("MuCapture/anDnLateResponse"+self.ResponseName) == None:
			self.unfoldIsMC = False
		else:
			self.unfoldIsMC = True

		self.unfoldresult = None

		self.nbiterations = 10
		# Range where the results are trusted
		self.minBin = 9
		self.maxBin = 17


	def InitTraining(self, thefile, defaulttitle):
		(self.training, self.trainingTitle) = Openfile(thefile, defaulttitle)

	
	def PrintMeasured(self):

		if self.logy:
			self.Canv.SetLogy()
		def DrawMeasuredFromResponse(thefile, thetitle, legend):
			print "MuCapture/anDnLateResponse"+self.ResponseName
			response = thefile.Get("MuCapture/anDnLateResponse"+self.ResponseName)
			recomeas = response.Hmeasured()
			if not thetitle == None:
				recomeas.SetTitle(thetitle+" measured spectrum")
			recomeas.Draw()
			Maxreco = recomeas.GetMaximum()

			print " ----> MuCapture/LateResponse"+self.ResponseName+"/MCTruthMomentumReco"
			truemeas = thefile.Get("MuCapture/LateResponse"+self.ResponseName+"/MCTruthMomentumReco")
			truemeas.SetLineColor(2)
			Maxtrue = truemeas.GetMaximum()

			MaxBin = max(Maxreco, Maxtrue)
			# recomeas.SetAxisRange(-(0.1*MaxBin),(1.1*MaxBin),'Y')
			recomeas.SetAxisRange(1e-6,(1.1*MaxBin),'Y')
			
			legend.AddEntry(recomeas, 'Reconstructed','L')
			legend.AddEntry(truemeas, 'Truth for fitted trks', 'L')
			truemeas.Draw("same")
		

		leg = TLegend(0.65,0.7,0.98,0.9)
		DrawMeasuredFromResponse(self.training, self.trainingTitle, leg)
		leg.Draw()
		self.Canv.Print(self.output)
		self.Canv.Clear()

		if self.unfoldIsMC:
			leg.Clear()
			DrawMeasuredFromResponse(self.unfold, self.unfoldTitle, leg)
			leg.Draw()
			self.Canv.Print(self.output)
			self.Canv.Clear()

		if self.logy:
			self.Canv.SetLogy(0)


	def PrintResponseMatrix(self):

		def DrawMeasuredFromResponse(thefile, thetitle):
			matrix = thefile.Get("MuCapture/LateResponse"+self.ResponseName+"/MCMeasVsTruthMomentum")
			if not thetitle == None:
				matrix.SetTitle(thetitle+" response matrix")
			matrix.Draw("colz")
		
		DrawMeasuredFromResponse(self.training, self.trainingTitle)
		self.Canv.Print(self.output)
		self.Canv.Clear()



	def UnfoldBayes(self):
		response = self.training.Get("MuCapture/anDnLateResponse"+self.ResponseName)

		# TODO: This is temporary until MeasuredMomentum is filled for both data and MC
		tmpresponse = self.unfold.Get("MuCapture/anDnLateResponse"+self.ResponseName)
		unfoldMeas = tmpresponse.Hmeasured()

		self.unfoldresult = RooUnfoldBayes(response, unfoldMeas, self.nbiterations);
		self.unfoldname = "Bayes"


	def PrintResultsMC(self):
		if self.unfoldresult == None:
			print "ERROR: No unfolding performed yet !!!"
			sys.exit(1)

		if self.logy:
			self.Canv.SetLogy()
		hReco= self.unfoldresult.Hreco()

		if self.unfoldIsMC:
			unfoldTruth = self.unfold.Get("MuCapture/LateResponse"+self.ResponseName+"/MCTruthMomentum")
			unfoldTruth.SetLineColor(2)
			leg = TLegend(0.65,0.7,0.98,0.9)

			MaxhReco = max(hReco.GetMaximum(), unfoldTruth.GetMaximum())
			MinhReco = 99999.
			for b in range(1,hReco.GetNbinsX()+1):
				Min = hReco.GetBinContent(b)
				if Min > 0.0:
					MinhReco = min(Min,MinhReco)
			unfoldTruth.Draw()
			hReco.Draw("same")
			leg.AddEntry(unfoldTruth, 'True input spectrum','L')
			leg.AddEntry(hReco, 'Unfolded spectrum', 'L')
			unfoldTruth.SetTitle("%s spectrum unfolded with %s spectrum using RooUnfold%s" % (self.unfoldTitle, self.trainingTitle, self.unfoldname))
			unfoldTruth.SetAxisRange((0.1*MinhReco),(1.1*MaxhReco),'Y')
			leg.Draw()
		else:
			hReco.Draw()
			hReco.SetTitle(self.trainingTitle+" spectrum unfolded using RooUnfold"+self.unfoldname)
		self.Canv.Update()
		self.Canv.Print(self.output)
		self.Canv.Clear()
		if self.logy:
			self.Canv.SetLogy(0)


	def PrintUnfoldDiff(self, Unfold1, Unfold2, Title1, Title2):
		if self.logy:
			self.Canv.SetLogy()

		hReco1= Unfold1.Hreco()
		hReco2= Unfold2.Hreco()

		hDiff = hReco1.Clone()
		hDiff.Add(hReco2, -1.)
		hDiff.Divide(hReco2)

		MinhReco = 999999.
		MaxhReco = -999999.
		for b in range(hDiff.GetNbinsX()):
			binCont = hDiff.SetBinError(b,0.0000001)
		for b in range(self.minBin,self.maxBin):
			binCont = hDiff.GetBinContent(b)
			MinhReco = min(MinhReco, binCont)
			MaxhReco = max(MaxhReco, binCont)
		if MinhReco < 0.0 and MaxhReco < 0.0:
			hDiff.SetAxisRange((1.1*MinhReco),(-0.1*MaxhReco),'Y')
		else:
			hDiff.SetAxisRange((1.1*MinhReco),(1.1*MaxhReco),'Y')
		hDiff.SetTitle("Relative difference %s minus %s for %s;Momentum [MeV/c]" % (Title1, Title2, self.unfoldTitle))
		hDiff.Draw()

		self.Canv.Update()
		self.Canv.Print(self.output)
		self.Canv.Clear()
		if self.logy:
			self.Canv.SetLogy(0)



	def MakeChi2Vs_Plot(self):
		self.RUChi2 = TGraph()
		self.MyChi2 = TGraph()

	def AddChi2Vs_Point(self, xVal):
		if self.unfoldresult == None:
			print "ERROR: No unfolding performed yet !!!"
			sys.exit(1)

		if not self.unfoldIsMC:
			print "WARNING: There is no truth information to compare to. No Chi2 point added."
			return

		unfoldTruth = self.unfold.Get("MuCapture/LateResponse"+self.ResponseName+"/MCTruthMomentum")
		chi2 = self.unfoldresult.Chi2(unfoldTruth)
		nbPt = self.RUChi2.GetN()
		self.RUChi2.SetPoint(nbPt, xVal, chi2)

		chi2 = 0
		hReco= self.unfoldresult.Hreco()
		for b in range(self.minBin,self.maxBin):
			yreco = hReco.GetBinContent(b)
			eyreco = hReco.GetBinError(b)
			ytrue = unfoldTruth.GetBinContent(b)
			chi2 = chi2 + ( (yreco-ytrue) * (yreco-ytrue) / (eyreco*eyreco))
		nbPt = self.MyChi2.GetN()
		self.MyChi2.SetPoint(nbPt, xVal, chi2)



	def DrawChi2Vs_Plot(self, xTitle):
		self.Canv.Clear()
		self.RUChi2.SetTitle("RooUnfold chi^{2} between unfolded and true %s spectra, using %s;%s;chi^{2}" % (self.unfoldTitle,self.trainingTitle,xTitle))
		self.RUChi2.SetMarkerStyle(5)
		self.RUChi2.SetMarkerColor(4)
		self.RUChi2.SetLineColor(4)
		self.RUChi2.Draw("APL")
		self.Canv.Update()
		self.Canv.Print(self.output)
		self.Canv.Clear()

		self.Canv.Clear()
		self.MyChi2.SetTitle("Manual chi^{2} between unfolded and true %s spectra, using %s;%s;chi^{2}" % (self.unfoldTitle,self.trainingTitle,xTitle))
		self.MyChi2.SetMarkerStyle(5)
		self.MyChi2.SetMarkerColor(4)
		self.MyChi2.SetLineColor(4)
		self.MyChi2.Draw("APL")
		self.Canv.Update()
		self.Canv.Print(self.output)
		self.Canv.Clear()

	def Close(self):
		self.Canv.Print(self.output+"]")



# if __name__ == "__main__":


# OptionParser.format_description = lambda self, formatter: "\n".join(map(textwrap.fill, self.description.split("\n"))) + "\n"
# parser = OptionParser(description=description)
parser = OptionParser()
parser.add_option("-t", "--training", dest="training", help="Root file containing the rooUnfold response.", default="")
parser.add_option("-u", "--unfold", dest="unfold", help="Root file containing the measured spectrum you want to unfold.")
parser.add_option("-o", "--output", dest="output", help="Output file name for the pdf file. Don't put the suffix '.pdf'.", default="output")
parser.add_option("", "--logy", dest="logy", help="Use SetLogY for the 1D plots.", action="store_true")
parser.add_option("-i", "--iterationcheck", dest="iterationcheck", help="Plot the unfolding results versus the number of iterations.", action="store_true")
parser.add_option("-m", "--multipletraining", dest="multipletraining", help="Instead of one training, loop over the training files given as arguments.", action="store_true")
parser.add_option("-d", "--differencetraining", dest="differencetraining", help="Perform the unfolding for two provided training files and plot the difference, first minus last file.", action="store_true")
parser.add_option("-r", "--responsename", dest="responsename", help="", default="")
(options, args) = parser.parse_args()

gROOT.SetBatch()
gStyle.SetOptStat(10)
gStyle.SetTitleW(1.0)

Unf = Unfolder(options)
# Study of unfolding results vs number of iterations
if options.iterationcheck:
	Unf.MakeChi2Vs_Plot()
	for it in range(1,16):
		Unf.nbiterations = it
		Unf.UnfoldBayes()
		Unf.unfoldname = "Bayes, %d iterations" % it
		Unf.PrintResultsMC()
		Unf.AddChi2Vs_Point("Number iterations")
	Unf.trainingTitle = Unf.trainingTitle+" spectrum"
	Unf.DrawChi2Vs_Plot()
elif options.multipletraining:
	Unf.MakeChi2Vs_Plot()
	for tr,TR in enumerate(args):
		Unf.InitTraining(TR, "training "+str(tr+1))
		Unf.UnfoldBayes()
		Unf.PrintResultsMC()
		Unf.AddChi2Vs_Point(tr+1)
	Unf.trainingTitle = str(len(args))+" training spectra"
	Unf.DrawChi2Vs_Plot("Training index")
elif options.differencetraining:
	if not len(args) == 2:
		print "ERROR: you did not provide two files in the command line."
		sys.exit()
	Unf.InitTraining(args[0], "1st training")
	Unf.UnfoldBayes()
	Unf.PrintResultsMC()
	UnfoldResults1 = Unf.unfoldresult
	TrainingTitle1 = Unf.trainingTitle

	Unf.InitTraining(args[1], "2nd training")
	Unf.UnfoldBayes()
	Unf.PrintResultsMC()
	UnfoldResults2 = Unf.unfoldresult
	TrainingTitle2 = Unf.trainingTitle

	Unf.PrintUnfoldDiff(UnfoldResults1, UnfoldResults2, TrainingTitle1, TrainingTitle2)

# Default output
else:
	Unf.PrintMeasured()
	Unf.PrintResponseMatrix()
	Unf.UnfoldBayes()
	Unf.PrintResultsMC()
Unf.Close()
