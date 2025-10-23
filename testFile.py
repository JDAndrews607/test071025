###Test data taken from IC1061, calibration curve from IC1056###


from time import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px     #for basic chromatogram plotting
import plotly.graph_objects as go   #for chromatogram curve fitting overlap
import csv
from random import random, randint, uniform, sample
import math
import scipy

eulersNumber = 2.718281828459045
pi = 3.141592653589793

#-------------------------------------------------------------#
#                 Calibration and Sample Data                 #
#-------------------------------------------------------------#

class SampleValues:
    def __init__(self, name, retentionTime = 0.0, lowerBound = 0.0, upperBound = 0.0):
        self.name = name
        #self.retentionTime = retentionTime
        #self.lowerBound = lowerBound
        #self.upperBound = upperBound
        print(f"{name} : retention time = {retentionTime} : counting between {lowerBound} and {upperBound}")
        #self.defaultCalibrationStandards()
        self.calculateCalibrationStandards()
        print(f"Calibration curve values: {self.calibration_Curve}")

    def calculateCalibrationStandards(self):
        for i in range(10):
            timeStamp, CD_signal = importChromatogram("../Calibration/"+str(i))
            corrected_CD = baselineCorrection_ALS(CD_signal)
            smoothed_signal = smoothing_SavitzkyGolay(corrected_CD, 27, 2)
            self.calibration_Curve.append(i)    #change to 


    def defaultCalibrationStandards(self):
        if self.name == "Li":
            self.retentionTime = 5.2
            self.calibration_Curve = [0.0066, 0.0621, 0.1232, 0.1844, 0.2459, 0.4919, 0.9841, 1.4764, 1.9675, 2.4596]
        elif self.name == "Na":
            self.retentionTime = 7.05
            self.calibration_Curve = [0.0378, 0.3558, 0.7055, 1.0564, 1.4085, 2.8176, 5.6377, 8.4575, 11.2706, 14.0896]
        elif self.name == "K":
            self.retentionTime = 9.0
            self.calibration_Curve = [0.1074, 1.0100, 2.0028, 2.9992, 3.9987, 7.9990, 16.0051, 24.0106, 31.9968, 39.9998]
        elif self.name == "Mg":
            self.retentionTime = 28.6
            self.calibration_Curve = [0.0134, 0.1258, 0.2495, 0.3736, 0.4980, 0.9963, 1.9934, 2.9905, 3.9852, 4.9820]
        elif self.name == "Ca":
            self.retentionTime = 42.2
            self.calibration_Curve = [0.0133, 0.1252, 0.2483, 0.3718, 0.4957, 0.9915, 1.9839, 2.9762, 3.9662, 4.9582]
        else:
            print(f"Component {self.name} is not supported")
            exit

    def getName(self):
        return self.name
    
    def get_RT(self):
        return self.retentionTime
    
    def get_lowerBound(self):
        return self.lowerBound

    def  get_upperBound(self):
        return self.upperBound
    
    def get_calibrationCurve(self):
        return self.calibration_Curve
    
    def setRetentionTime(self, newTime):
        self.retentionTime = newTime


    def findPeakBoundaries_TLM(self, timestamp, CD_signal): #tangent line method
        pass

    #For compound
    name = " "
    retentionTime = 0.0
    lowerBound = 0.0
    upperBound = 0.0

    #For calibration curve:
    calibration_Curve = []




#-------------------------------------------------------------#
#                        Test Data                            #
#-------------------------------------------------------------#

def generate_IC_NoisyData(datapoints = 10000, noise = 0.05):
    time = np.linspace(0, 60, datapoints)

    def genGaussianPeak(x, mu, sigma, amplitude):
        curve = amplitude * eulersNumber**(-0.5 * ((x - mu) / sigma)**2)
        return curve

    signal = (
        genGaussianPeak(time, 7, 0.5, 0.25) +
        genGaussianPeak(time, 9, 0.75, 0.2) +
        genGaussianPeak(time, 11, 0.5, 0.3) +
        genGaussianPeak(time, 24, 0.5, 0.2) +
        genGaussianPeak(time, 42, 0.20, 0.2)
    )

    baselineDrift = (
        0.001 * np.linspace(0, 1, len(signal)) +
        0.002 * np.linspace(-1, 1, len(signal))**2 +
        0.005 * (1 - np.linspace(0, 5, len(signal)))
    )
    signal += baselineDrift

    noise = np.random.normal(0, noise, datapoints)
    noiseySignal = signal + noise

    return time.tolist(), noiseySignal.tolist()

#-------------------------------------------------------------#
#                    Chromatogram Tools                       #
#-------------------------------------------------------------#

def importChromatogram(fileLocation):
    timestamp = []
    CD_signal = []

    with open(fileLocation, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for row in reader:
            #print("Time:", row[0], ", CD Signal:", row[1])
            timestamp.append(float(row[0].lstrip('\ufeff')))
            CD_signal.append(float(row[1].lstrip('\ufeff')))
        
    return timestamp, CD_signal


def plotChromatogram(timestamp, CD_signal):
    dataframe = pd.DataFrame({'Retention Time': timestamp, 'CD': CD_signal})
    fig = px.line(dataframe, x='Retention Time',y='CD', title='Chromatogram')
    fig.show()


def plotChromatogram_EMG_overlap(original_timestamp, original_CD_Signal, EMG_timestamp, EMG_CD_Signal):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=original_timestamp, y=original_CD_Signal, mode='lines', name='Original Chromatogram'))
    fig.add_trace(go.Scatter(x=EMG_timestamp, y=EMG_CD_Signal, mode='lines', name='EMG Overlap'))
    fig.update_layout(
        xaxis_title = 'Retention Time (overlap test)',
        yaxis_title = 'CD Signal'
    )
    fig.show()

#-------------------------------------------------------------#
#                      Curve Smoothing                        #
#-------------------------------------------------------------#

def smoothing_SavitzkyGolay(CD_signal, windowSize, polynomialOrder):
    if windowSize % 2 == 0 or windowSize <= polynomialOrder:
        raise ValueError("Window size must be odd and greater than the order of polynomial")

    halfWindow = windowSize // 2    #double // to do int division instead of float
    size = polynomialOrder + 1

    #Building matrix A
    x = [i - halfWindow for i in range(windowSize)] #provided centered window positions
    A = [] #Vandermonde matrix
    for i in range(windowSize):
        row = []
        for j in range(size):
            row.append(x[i]**j)
        A.append(row)

    #Computing transpose A * A (A^T*A)
    A_Transposed_A = []
    for i in range(size):
        row = []
        for j in range(size):
            s = 0.0
            for k in range(windowSize):
                s+= A[k][i] * A[k][j]
            row.append(s)
        A_Transposed_A.append(row)
    
    #Computing Transpose A (A^T)
    A_Transposed = []
    for i in range(size):
        row = []
        for j in range(windowSize):
            row.append(A[j][i])
        A_Transposed.append(row)

    #Compute pseudoinverse then eliminate w/ Gauss Jordan
    aug_matrix = []
    for i in range(size):
        aug_matrix.append(A_Transposed_A[i] + [1.0 if i==j else 0.0 for j in range(size)])

    for i in range(size):
        pivot = aug_matrix[i][i]
        if abs(pivot) < 1e-12:
            raise ValueError("Singular Matrix")
        for j in range(2 * size):
            aug_matrix[i][j] /= pivot
        for k in range(size):
            if k!=i:
                factor = aug_matrix[k][i]
                for j in range(2 * size):
                    aug_matrix[k][j] -= factor * aug_matrix[i][j]

    ATA_Inverse = []
    for i in range(size):
        ATA_Inverse.append(aug_matrix[i][size:])

    #Calculate pseudoinverse
    psuedoInverse = []
    for i in range(size):
        row = []
        for j in range(windowSize):
            s = 0.0
            for k in range(size):
                s += ATA_Inverse[i][k] * A_Transposed[k][j]
            row.append(s)
        psuedoInverse.append(row)

    #Obtain smoothing coefficients
    sm_Coeff = psuedoInverse[0]

    #Apply convolutional with edge padding
    smoothed_signal = []
    padding = [CD_signal[0]]*halfWindow + CD_signal + [CD_signal[-1]] * halfWindow
    for i in range(len(CD_signal)):
        window = padding[i:i + windowSize]
        smoothedValue = 0.0
        for j in range(windowSize):
            smoothedValue += sm_Coeff[j]*window[j]
        smoothed_signal.append(smoothedValue)

    return smoothed_signal


#use Thomas algorithm for tridiagonal system
def baselineCorrection_ALS(CD_signal):   #asymetric least squares
    length = len(CD_signal)
    weight = [1.0] * length
    baseline = [0.0] * length
    _lambda = 1e5   #smoothing parameter
    assymParam_p = 0.01     #assymetry parameter
    old_Baseline = None   ##added to allow for early stop if negligable baseline change between iterations

    #components for tridiagonal matrix
    a = [_lambda] * (length - 1)
    b_n = [1.0 + 2 * _lambda] * length
    c = [_lambda] * (length - 1)

    count = 10
    for _ in range(count):
        b = []
        d = []
        for i in range(length):
            b.append(b_n[i] + weight[i] - 1)
            d.append(weight[i] * CD_signal[i])

        #Thomas Algorithm
        for i in range(1, length):
            m = a[i-1] / b[i-1]
            b[i] -= m * c[i-1]
            d[i] -= m * d[i-1]
        baseline = [0.0] * length
        baseline[-1] = d[-1] / b[-1]
        for i in reversed(range(length - 1)):
            baseline[i] = (d[i] - c[i] * baseline[i+1]) / b[i]

        #updating weights
        for i in range(length):
            if CD_signal[i] > baseline[i]:
                weight[i] = assymParam_p
            else:
                weight[i] = 1 - assymParam_p

    corrected_CD = []
    for i in range(length):
        corrected_CD.append(CD_signal[i] - baseline[i])
    return corrected_CD

def baselineCorrection_ShirleyTougaard():
    pass
    

#-------------------------------------------------------------#
#                Calibration Curve Determination              #
#-------------------------------------------------------------#


def compareStandardToTrue(standardVal, trueVal):
    #Dictonaries to hold which standards work best for each triple and pair of standards
    standardStoreThreeVal = {}
    standardStoreTwoVal = {}
    finalRanges = {}
    #Standard ranges will be compared to find best combination of standards and integration formulas for best reults
    for i in range(3,11):
        print(f"Testing standards {i-2}, {i-1}, and {i}")
        for j in range(i-2, i+1):
            print(f"Comparing standard {j}")
            key = f"{i-2}{i-1}{i}_{j}"
            methods = {
                'linNoOffset' : (integrate_LinearNoOffset(standardVal[j-1]) - integrate_LinearNoOffset(trueVal[j-1]))**2,
                '1OverArea' : (integrate_1OverArea(standardVal[j-1]) - integrate_1OverArea(trueVal[j-1]))**2,
                '1OverAreaSquared' : (integrate_1OverLogAreaSquared(standardVal[j-1]) - integrate_1OverLogAreaSquared(trueVal[j-1]))**2,
                '1OverResponse' : (integrate_1OverResponse(standardVal[j-1]) - integrate_1OverResponse(trueVal[j-1]))**2,
                '1OverResponseSquared' : (integrate_1OverResponseSquared(standardVal[j-1]) - integrate_1OverResponseSquared(trueVal[j-1]))**2,
                '1OverLogArea' : (integrate_1OverLogArea(standardVal[j-1]) - integrate_1OverLogArea(trueVal[j-1]))**2,
                '1OverLogAreaSquared' : (integrate_1OverLogAreaSquared(standardVal[j-1]) - integrate_1OverLogAreaSquared(trueVal[j-1]))**2,
                '1OverLogResponse' : (integrate_1OverLogResponse(standardVal[j-1]) - integrate_1OverLogResponse(trueVal[j-1]))**2,
                '1OverLogResponseSquared' : (integrate_1OverLogResponseSquared(standardVal[j-1]) - integrate_1OverLogResponseSquared(trueVal[j-1]))**2
            }
            best_method = min(methods, key = methods.get)
            best_value = methods[best_method]
            standardStoreThreeVal[key] = {
                'method' : best_method,
                'value' : best_value
            } 
    for key, result in standardStoreThreeVal.items():
        print(f"Standards {key}: Best Method = {result['method']}: Mean Square Difference = {result}")
    
    #Gives which pair of standards provide best fit for given range
        #loop over standards and create pairwise sets for later comparison
        #Done by comparing pairs of standards and finding which values have lowest mean square distance in values
    for i in range(2,11):       #DO NOT CHANGE. WILL BREAK PROGRAM!
        for j in range(0,2):
            pairwise_key = f"stdRange_{i-1}{i}{i+1}_{i-1+j}_{i+j}"
            if i > 8:
                standardStoreTwoVal.pop("stdRange_91011_9_10", None)    #removes extraneous standard range if present
                standardStoreTwoVal.pop("stdRange_91011_10_11", None)   #removes extraneous standard range if present   (which should be eventually after i = 7)

            if i+1 < 11:    #prevents accessing extraneous standard ranges from previous loop
                key_i = f"{i-1}{i}{i+1}_{i}"
                key_iMinus1 = f"{i-1}{i}{i+1}_{i-1}"

                val_i = standardStoreThreeVal[key_i]['value']
                val_iMinus1 = standardStoreThreeVal[key_iMinus1]['value']

                best_key = key_i if val_i < val_iMinus1 else key_iMinus1
                best_method = standardStoreThreeVal[best_key]['method']
                best_value = standardStoreThreeVal[best_key]['value']

            standardStoreTwoVal[pairwise_key] = {
                'method' : best_method,
                'value' : best_value
            }
    standardStoreTwoVal.pop("stdRange_91011_10_11", None)   #removes standard range that shouldn't be present. Not sure why it gets made...

    #makes final dictionary list of standard ranges/integration formulas
    finalRanges = {
        "stdRange_1_2" : {
            'method' : standardStoreTwoVal['stdRange_123_1_2']['method'],
            'value' : 100 - standardStoreTwoVal['stdRange_123_1_2']['value'] *100
        },
        
        "stdRange_2_3" : {
            'method' : standardStoreTwoVal['stdRange_123_2_3']['method'] if standardStoreTwoVal['stdRange_123_2_3']['value'] < standardStoreTwoVal['stdRange_234_2_3']['value'] else standardStoreTwoVal['stdRange_234_2_3']['method'],
            'value' : 100 - (standardStoreTwoVal['stdRange_123_2_3']['value'] if standardStoreTwoVal['stdRange_123_2_3']['value'] < standardStoreTwoVal['stdRange_234_2_3']['value'] else standardStoreTwoVal['stdRange_234_2_3']['value']) * 100
        },

        "stdRange_3_4" : {
            'method' : standardStoreTwoVal['stdRange_234_3_4']['method'] if standardStoreTwoVal['stdRange_234_3_4']['value'] < standardStoreTwoVal['stdRange_345_3_4']['value'] else standardStoreTwoVal['stdRange_345_3_4']['method'],
            'value' : 100 - (standardStoreTwoVal['stdRange_234_3_4']['value'] if standardStoreTwoVal['stdRange_234_3_4']['value'] < standardStoreTwoVal['stdRange_345_3_4']['value'] else standardStoreTwoVal['stdRange_345_3_4']['value']) * 100
        },

        "stdRange_4_5" : {
            'method' : standardStoreTwoVal['stdRange_345_4_5']['method'] if standardStoreTwoVal['stdRange_345_4_5']['value'] < standardStoreTwoVal['stdRange_456_4_5']['value'] else standardStoreTwoVal['stdRange_456_4_5']['method'],
            'value' : 100 - (standardStoreTwoVal['stdRange_345_4_5']['value'] if standardStoreTwoVal['stdRange_345_4_5']['value'] < standardStoreTwoVal['stdRange_456_4_5']['value'] else standardStoreTwoVal['stdRange_456_4_5']['value']) * 100
        },

        "stdRange_5_6" : {
            'method' : standardStoreTwoVal['stdRange_456_5_6']['method'] if standardStoreTwoVal['stdRange_456_5_6']['value'] < standardStoreTwoVal['stdRange_567_5_6']['value'] else standardStoreTwoVal['stdRange_567_5_6']['method'],
            'value' : 100 - (standardStoreTwoVal['stdRange_456_5_6']['value'] if standardStoreTwoVal['stdRange_456_5_6']['value'] < standardStoreTwoVal['stdRange_567_5_6']['value'] else standardStoreTwoVal['stdRange_567_5_6']['value']) * 100
        },

        "stdRange_6_7" : {
            'method' : standardStoreTwoVal['stdRange_567_6_7']['method'] if standardStoreTwoVal['stdRange_567_6_7']['value'] < standardStoreTwoVal['stdRange_678_6_7']['value'] else standardStoreTwoVal['stdRange_678_6_7']['method'],
            'value' : 100 - (standardStoreTwoVal['stdRange_567_6_7']['value'] if standardStoreTwoVal['stdRange_567_6_7']['value'] < standardStoreTwoVal['stdRange_678_6_7']['value'] else standardStoreTwoVal['stdRange_678_6_7']['value']) * 100
        },

        "stdRange_7_8" : {
            'method' : standardStoreTwoVal['stdRange_678_7_8']['method'] if standardStoreTwoVal['stdRange_678_7_8']['value'] < standardStoreTwoVal['stdRange_789_7_8']['value'] else standardStoreTwoVal['stdRange_789_7_8']['method'],
            'value' : 100 - (standardStoreTwoVal['stdRange_678_7_8']['value'] if standardStoreTwoVal['stdRange_678_7_8']['value'] < standardStoreTwoVal['stdRange_789_7_8']['value'] else standardStoreTwoVal['stdRange_789_7_8']['value']) * 100
        },

        "stdRange_8_9" : {
            'method' : standardStoreTwoVal['stdRange_789_8_9']['method'] if standardStoreTwoVal['stdRange_789_8_9']['value'] < standardStoreTwoVal['stdRange_8910_8_9']['value'] else standardStoreTwoVal['stdRange_8910_8_9']['method'],
            'value' : 100 - (standardStoreTwoVal['stdRange_789_8_9']['value'] if standardStoreTwoVal['stdRange_789_8_9']['value'] < standardStoreTwoVal['stdRange_8910_8_9']['value'] else standardStoreTwoVal['stdRange_8910_8_9']['value']) * 100
        },

        "stdRange_9_10" : {
            'method' : standardStoreTwoVal['stdRange_8910_9_10']['method'],
            'value' : 100 - standardStoreTwoVal['stdRange_8910_9_10']['value'] * 100
        }
    }

    return finalRanges  #this will be a dictionary
    

#-------------------------------------------------------------#
#                    Peak Area Calculation                    #
#-------------------------------------------------------------#


def peakFind_EMG(timestamp, CD_Signal, targetTime, _lambda = 0.01, peakDeviation = 1):     #Using exponentially modified Gaussian equation
    #lambda is the rate of the exponential component
    #set peakDeviation to how far away from the targetTime you allow the program to check for the peak. Setting to 5 will be five minutes +/- targetTime

    #find peak close to the estimated peak for your ion
    maxPeak = 0.0
    peakLocation = 0.0
    peakIterator = 0
    for i in range(len(timestamp)):
        if CD_Signal[i] > maxPeak and (timestamp[i] > targetTime - peakDeviation and timestamp[i] < targetTime + peakDeviation):
            maxPeak = CD_Signal[i]
            peakLocation = timestamp[i]
            peakIterator = i
    print(f"Peak located at {peakLocation} with signal {maxPeak}")

    lowerBound = 0.0
    lowerIterator = 0
    upperBound = 0.0
    upperIterator = 0
    #find boundaries of this peak
        #lower bound:
    peakScalar = 0.01 #determines what percent of peak height is used to end peak
    for i in range(peakIterator, 0, -1):
        if CD_Signal[i] <= peakScalar * maxPeak:  
            lowerBound = timestamp[i+1]
            lowerIterator = i+1
            break
        #upper bound
    for i in range(peakIterator, len(timestamp)-1):
        if CD_Signal[i] <= peakScalar * maxPeak:
            upperBound = timestamp[i-1]
            upperIterator = i-1
            break
    print(f"Peak bound between {lowerBound}, {CD_Signal[lowerIterator]} and {upperBound}, {CD_Signal[upperIterator]}")

    #calculate mean and standard deviation
    mu = 0.0    #center of the peak, computed as central mass of peak
    weight = 0.0
    total = 0.0
    for i in range(lowerIterator, upperIterator):
        weight += timestamp[i] * CD_Signal[i]
        total += CD_Signal[i]
    mu = weight / total
    print(f"Mean = {mu}")

    sigma = 0.0 #standard deviation, using weighted variance
    variance = 0.0
    for i in range(lowerIterator, upperIterator):
        variance += CD_Signal[i] * (timestamp[i] - mu)**2
    sigma = (variance / total)**(1/2)   #reusing total variable from mu calculation
    print(f"Standard Deviation = {sigma}")

    #save peak and boundary information
    peak_timestamp = timestamp[lowerIterator:upperIterator+1]
    peak_CD_Signal = CD_Signal[lowerIterator:upperIterator+1]

    #approximate area with Riemann sums:
    area = 0.0
    for i in range(len(peak_CD_Signal) - 1):
        dt = peak_timestamp[i+1] - peak_timestamp[i]
        area += peak_CD_Signal[i] * dt

    #fit to EMG peak and optimize values
    area, mu, sigma, _lambda = fit_EMG(peak_timestamp, peak_CD_Signal, area, mu, sigma, _lambda)

    #calculate EMG values for each signal point
    updated_CD_Signal = []
    for i in peak_CD_Signal:
        updated_CD_Signal.append(exponentialModifiedGaussian_EMG(i, mu, sigma, _lambda, area))

    plotChromatogram_EMG_overlap(timestamp, CD_Signal, peak_timestamp, updated_CD_Signal)
    return area

#sum of squared errors. Feed into curve fit to 
def errorFunction_EMG(timestamp, CD_Signal, area, mu, sigma, _lambda):
    error = 0.0
    for i in range(len(CD_Signal)):
        modelVal = exponentialModifiedGaussian_EMG(CD_Signal[i],  mu, sigma, _lambda, area)
        error += (CD_Signal[i] - modelVal)**2
    return error


#Two step fitting method to determine best parameters. Using differential evolution followed by gradient descent
def fit_EMG(timestamp, CD_Signal, init_area, init_mu, init_sigma, init_lambda, steps = 300, step_size = 0.01, learning_Rate = 0.001):
    #Global method to narrow down correct parameters
    #Using differential evolution method

        #setting arrays for use in differential convolution
    np.random.seed(42)
    diffEvo_Params = []
    bounds = [  #bounded values for each parameter
        (init_mu - 5*init_sigma, init_mu + 5*init_sigma),     #for mu
        (0.1 * init_sigma, 5 * init_sigma),   #for sigma
        (0.001, 1.0),   #for lambda
        (0.1 * init_area, 10 * init_area)   #for area
    ]
    
    population_size = 15
    dimension = len(bounds)
    
    population = []
    for i in range(population_size):
        individual = []
        for j in range(dimension):
            individual.append(uniform(*bounds[j]))  #using random.uniform
        population.append(individual)

    scores = []
    for i in population:
        mu, sigma, _lambda, area = i
        val = 0
        for j in range(len(CD_Signal)):
            y_pred = exponentialModifiedGaussian_EMG(CD_Signal[j], mu, sigma, _lambda, area)
            val += (y_pred - CD_Signal[j])**2  #using weighted averages
        val /= len(CD_Signal)
        scores.append(val)

    #Setting up sequence for differential evolution
    for step in range(0, steps):    #1000 generations
        for i in range(population_size):
            index = list(range(population_size))
            index.remove(i)
            r1, r2, r3 = sample(index, 3)    #randomly chooses three sample values from the population. #using random.sample

            #Mutation
            mutant = []
            for j in range(dimension):
                val = population[r1][j] + step_size * (population[r2][j] - population[r3][j])
                mutant.append(val)
                #clamp the mutant list
            for j in range(dimension):
                mutant[j] = max(bounds[j][0], min(mutant[j], bounds[j][1]))

            #Crossover
            trial = []
            crossoverProb = 0.7
            for j in range(dimension):
                if random() < crossoverProb or j == randint(0, dimension - 1):
                    trial.append(mutant[j])
                else:
                    trial.append(population[i][j])

            #evaluate crossover trial
            mu = trial[0]
            sigma = trial[1]
            _lambda = trial[2]
            area = trial[3]
            trialError = 0
            for j in range(len(CD_Signal)):
                y_pred = exponentialModifiedGaussian_EMG(CD_Signal[j], mu, sigma, _lambda, area)
                trialError += (y_pred - CD_Signal[j])**2
            trialError /= len(CD_Signal)

            #Selection
            if trialError < scores[i]:
                population[i] = trial
                scores[i] = trialError

        bestIndex = scores.index(min(scores))
        diffEvo_Params = population[bestIndex]    #obtains the most accurate parameters from differential evolution [mu, sigma, lambda, area]
        

    #return diffEvo_Params

            #####End of Differential Evolution#####
            #Moving to Gradient Descent

    #local method to find final structure
    #using gradiant descent, starting parameters taken from differential evolution
    gd_mu = diffEvo_Params[0]
    gd_sigma = diffEvo_Params[1]
    gd_lambda = diffEvo_Params[2]
    gd_area = diffEvo_Params[3]
    count = 0
    for step in range(steps):
        delta = max(1e-5, 0.01 * abs(gd_mu))    #adaptively scale the learning parameter with mu
        new_mu = (errorFunction_EMG(timestamp, CD_Signal, gd_area, gd_mu+delta, gd_sigma, gd_lambda) - errorFunction_EMG(timestamp, CD_Signal, gd_area, gd_mu-delta, gd_sigma, gd_lambda)) / (2 * delta)
        gd_mu -= new_mu * learning_Rate
        gd_mu = max(gd_mu, 1e-6)
        new_sigma = (errorFunction_EMG(timestamp, CD_Signal, gd_area, gd_mu, gd_sigma+delta, gd_lambda) - errorFunction_EMG(timestamp, CD_Signal, gd_area, gd_mu, gd_sigma-delta, gd_lambda)) / (2 * delta)
        gd_sigma -= new_sigma * learning_Rate
        gd_sigma = max(gd_sigma, 1e-6)
        new_lambda = (errorFunction_EMG(timestamp, CD_Signal, gd_area, gd_mu, gd_sigma, gd_lambda+delta) - errorFunction_EMG(timestamp, CD_Signal, gd_area, gd_mu, gd_sigma, gd_lambda-delta)) / (2 * delta)
        gd_lambda -= new_lambda * learning_Rate
        gd_lambda = max(gd_lambda, 1e-6)
        new_area = (errorFunction_EMG(timestamp, CD_Signal, gd_area+delta, gd_mu, gd_sigma, gd_lambda) - errorFunction_EMG(timestamp, CD_Signal, gd_area-delta, gd_mu, gd_sigma, gd_lambda)) / (2 * delta)
        gd_area -= new_area * learning_Rate
        gd_area = max(gd_area, 1e-6)

        #convergence check
        if abs(new_mu) < 1e-6 and abs(new_sigma) < 1e-6 and abs(new_lambda) < 1e-6 and abs(new_area) < 1e-6:
            print("Converged Early")
            break

        count+=1
        if count%100 == 0:
            print(f"Final fit: mu={gd_mu}, sigma={gd_sigma}, lambda={gd_lambda}, area={gd_area}")
    
    return gd_area, gd_mu, gd_sigma, gd_lambda
   

#actual exponential modified gaussian equation
def exponentialModifiedGaussian_EMG(x, mu, sigma, _lambda, area = 1):
    EMG_Front = _lambda/2 * eulersNumber**(_lambda/2 * (2*mu + _lambda*sigma**2 - 2*x))
    #erf = 0     #error function. Obtained from series shown below
    #erfc = 0    #1-erf(erfc_term)
    erfc_term = (mu + _lambda*sigma**2 - x) / (2**(1/2) * sigma)
    #Use Power Series Expansion of the Error Function
    #count = 0   #counter for erf power series
    #finding erf(erfc_term) with a power series expansion
    #while count < 10: #takes series to ten terms:
    #    erf += ( (-1)**count * erfc_term**(2*count+1) ) / ( factorial(count) * (2*count+1) )
    #    count += 1
    #erfc = 1 - 2 / pi**(1/2) * erf
    #EMG_final = EMG_Front * erfc_Abramowitz__Stegun(erfc_term)      #straight line, slight downward trend left to right
    #EMG_final = EMG_Front * erfc        #definitely don't use
    #EMG_final = area * EMG_Front * math.erfc(erfc_term)        #only works on scalars, not arrays
    EMG_final = area * EMG_Front * scipy.special.erfc(erfc_term)   #to control for broadcast issue. If removed, get rid of include scipy in program header

    return EMG_final

def erfc_Abramowitz__Stegun(val):   #good for both large and small signals
    t = 1 / (1 + 0.5 * abs(val))
    tau = t * eulersNumber**(-val * val - 1.26551223 + val * (
        1.00002368 + t * (
        0.37409196 + t * (
        0.09678418 + t * (
        -0.18628806 + t * (
        0.27886807 + t * (
        -1.13520398 + t * (
        1.48851587 + t * (
        -0.82215223 + t * 0.17087277
        )))))))))

    if val >= 0:
        return tau
    else:
        return 2 - tau

#find factorial of provided value
def factorial(factor):
    val = 1
    for i in range(2, factor + 1):
        val *= i
    return val




#-------------------------------------------------------------#
#                       Peak Integration                      #
#-------------------------------------------------------------#

    #Calculate based off of calibration standards. This will 
    #create a curve to compare your sample's signal area to
    #to determine concentration

def integrate_LinearNoOffset(val):
    return np.random.uniform(0, 1)

def integrate_1OverArea(val):
    return np.random.uniform(0, 1)

def integrate_1OverAreaSquared(val):
    return np.random.uniform(0, 1)

def integrate_1OverResponse(val):
    return np.random.uniform(0, 1)

def integrate_1OverResponseSquared(val):
    return np.random.uniform(0, 1)

def integrate_1OverLogArea(val):
    return np.random.uniform(0, 1)

def integrate_1OverLogAreaSquared(val):
    return np.random.uniform(0, 1)

def integrate_1OverLogResponse(val):
    return np.random.uniform(0, 1)

def integrate_1OverLogResponseSquared(val):
    return np.random.uniform(0, 1)




###MAIN###

def importTest():
    timestamp, CD_signal = importChromatogram(f"Calibration/{1}.csv")
    plotChromatogram(timestamp, CD_signal)
    corrected_CD = baselineCorrection_ALS(CD_signal)
    plotChromatogram(timestamp, corrected_CD)
    smoothed_signal = smoothing_SavitzkyGolay(corrected_CD, 27, 2)
    plotChromatogram(timestamp, smoothed_signal)
    peakFind_EMG(timestamp, smoothed_signal, 5.0)


def randomGenTest():
    timestamp, CD_signal = generate_IC_NoisyData(10000, 0.05)
    plotChromatogram(timestamp, CD_signal)
    corrected_CD = baselineCorrection_ALS(CD_signal)
    plotChromatogram(timestamp, corrected_CD)
    counter = 3
    countMax = 150
    while counter <= countMax:
        print("Window size="+str(counter))
        smoothed_signal = smoothing_SavitzkyGolay(corrected_CD, counter, 2)
        plotChromatogram(timestamp, smoothed_signal)
        counter += int(countMax / 5)
        if counter % 2 == 0:
            counter += 1

def curveFitTest():
    timestamp, CD_signal = importChromatogram(f"Calibration/{1}.csv")
    plotChromatogram(timestamp, CD_signal)
    corrected_CD = baselineCorrection_ALS(CD_signal)
    plotChromatogram(timestamp, corrected_CD)
    smoothed_signal = smoothing_SavitzkyGolay(corrected_CD, 27, 2)
    plotChromatogram(timestamp, smoothed_signal)
    peakFind_EMG(timestamp, smoothed_signal, 5)


curveFitTest()

#testing
#Li = SampleValues("Li")
#randomGenTest()
#importTest()
#curveFitTest()

#compx = [0.2, 0.5, 1.1, 1.4, 1.6, 1.3, 1.0, 0.8, 0.5, 0.1]
#compy = [0.22, 0.51, 1.01, 1.1, 1.45, 1.1, 1.0, 0.2, 0.45, 0.11]
#standardChoices = compareStandardToTrue(compx, compy)
#for key, val in standardChoices.items():
#    print(f"{key} : {val['method']}, Certainty = {val['value']}")
