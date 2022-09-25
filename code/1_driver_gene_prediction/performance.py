import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt

def cal_ROAUC(X, Y, path=False):
    fpr, tpr, thresholds = metrics.roc_curve(Y, X, pos_label=1)
    AUC = metrics.auc(fpr, tpr)
    if path:
        f1 = plt.plot(fpr, tpr)
        f1.savefig(path)
    return AUC

def cal_PRAUC(X,Y, path=False):
    prec, reca, thresholds = metrics.precision_recall_curve(Y, X, pos_label=1)
    AUC = metrics.auc(reca, prec)
    if path:
        f1 = plt.plot(reca, prec)
        f1.savefig(path)
    return AUC

def get_TPFNTNFP_from_strip(actual_list, predicted_list):
	TP, FN, TN, FP = 0, 0, 0, 0
	for i in range(len(actual_list)):
		if actual_list[i] == 1:
			if predicted_list[i] == 1: TP += 1
			elif predicted_list[i] == 0: FN += 1
		elif actual_list[i] == 0:
			if predicted_list[i] == 1: FP += 1
			elif predicted_list[i] == 0: TN += 1
	return [TP, FN, TN, FP]

def sensitivity(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return TP/float(TP+FN)
	except: return np.nan

def specificity(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return TN/float(TN+FP)
	except: return np.nan

def precision(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return TP/float(TP+FP)
	except: return np.nan

def recall(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return TN/float(TN+FN)
	except: return np.nan

def accuracy(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return (TP+TN)/float(TP+FN+TN+FP)
	except: return np.nan

def balanced_accuracy(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return (sensitivity(TPFNTNFP) + specificity(TPFNTNFP))/2
	except: return np.nan

def F1_score(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return 2*TP/float(2*TP+FP+FN)
	except: return np.nan

def MCC(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return float(TP*TN-FP*FN)/(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))**0.5)
	except: return np.nan

def fpr(TPFNTNFP, cri=['num', 'list'][0]):
	TP, FN, TN, FP = TPFNTNFP
	if cri == 'list':
		TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
	if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
	try: return float(FP)/(FP+TN)
	except: return np.nan

def get_TPFNTNFP_num(X, Y, threshold, reverse=False):
	TP, FN, TN, FP = 0, 0, 0, 0
	if not reverse:
		for i in range(len(Y)):
			if Y[i] == 1:
				if X[i] > threshold: TP += 1
				else: FN += 1
			else:
				if X[i] > threshold: FP += 1
				else: TN += 1
	else:
		for i in range(len(Y)):
			if Y[i] == 1:
				if X[i] < threshold: TP += 1
				else: FN += 1
			else:
				if X[i] < threshold: FP += 1
				else: TN += 1
	return [TP, FN, TN, FP] 

def find_threshold(X,Y, reverse=False):
	value_list = list(set(X))
	value_list = sorted(value_list, reverse=reverse)
	F1_index, MCC_index, Acc_index, bAcc_index = [0, 0], [0, 0], [0, 0], [0, 0]
	for value in value_list:
		TPFNTNFP = get_TPFNTNFP_num(X,Y,value,reverse=reverse)
		Acc_value = accuracy(TPFNTNFP)
		bAcc_value = balanced_accuracy(TPFNTNFP)
		MCC_value = MCC(TPFNTNFP)
		F1_value = F1_score(TPFNTNFP)
		if F1_index[0] < F1_value: F1_index = [F1_value, value]
		if MCC_index[0] < MCC_value: MCC_index = [MCC_value, value]
		if Acc_index[0] < Acc_value: Acc_index = [Acc_value, value]
		if bAcc_index[0] < bAcc_value: bAcc_index = [bAcc_value, value]
	return F1_index, MCC_index, Acc_index, bAcc_index

def find_threshold_fpr(X,Y, reverse=False, cri_fpr=0.2):
	value_list = list(set(X))
	value_list = sorted(value_list, reverse=reverse)
	fpr_index = [1, 0]
	for value in value_list:
		TPFNTNFP = get_TPFNTNFP_num(X,Y,value,reverse)
		tmp_fpr = fpr(TPFNTNFP)
		if abs(fpr_index[0]-cri_fpr) > abs(tmp_fpr-cri_fpr): fpr_index = [tmp_fpr, value]
	return fpr_index

def find_threshold_tpr(X,Y, reverse=False, cri_tpr=0.8):
	value_list = list(set(X))
	value_list = sorted(value_list, reverse=reverse)
	tpr_index = [1, 0]
	for value in value_list:
		TPFNTNFP = get_TPFNTNFP_num(X,Y,value,reverse)
		tmp_tpr = sensitivity(TPFNTNFP)
		if abs(tpr_index[0]-cri_tpr) > abs(tmp_tpr-cri_tpr): tpr_index = [tmp_tpr, value]
	return tpr_index

def get_TPFNTNFP_num2(X1, X2, Y, threshold1, threshold2, reverse2=False):
	TP, FN, TN, FP = 0, 0, 0, 0
	if not reverse2:
		for i in range(len(Y)):
			if Y[i] == 1:
				if (X1[i] > threshold1 or X2[i] > threshold2): TP += 1
				else: FN += 1
			else:
				if (X1[i] > threshold1 or X2[i] > threshold2): FP += 1
				else: TN += 1
	else:
		for i in range(len(Y)):
			if Y[i] == 1:
				if (X1[i] > threshold1 or X2[i] < threshold2): TP += 1
				else: FN += 1
			else:
				if (X1[i] > threshold1 or X2[i] < threshold2): FP += 1
				else: TN += 1
	return [TP, FN, TN, FP] 

