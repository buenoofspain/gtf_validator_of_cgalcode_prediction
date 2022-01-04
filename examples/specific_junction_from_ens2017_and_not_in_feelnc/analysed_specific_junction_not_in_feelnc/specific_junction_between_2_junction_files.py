#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:53:31 2019

@author: niguilla
"""

dict_pred = dict()
dict_ref = dict()

list_pred = list()
list_pred_tr = list()
list_ref = list()
list_ref_tr = list()

nb_junction_pred = 0
nb_junction_ref = 0
nb_test = 0

list_83 = ['CGACAFTENST00000297620JUNC00002', 'CGACAFTENST00000542447JUNC00006', 'CGACAFTENST00000536727JUNC00007', 'CGACAFTENST00000309703JUNC00008', 'CGACAFTENST00000380138JUNC00012', 'CGACAFTENST00000457054JUNC00017', 'CGACAFTENST00000488267JUNC00020', 'CGACAFTENST00000542836JUNC00028', 'CGACAFTENSMUST00000100772JUNC00037', 'CGACAFTENST00000380036JUNC00041', 'CGACAFTENST00000404234JUNC00049', 'CGACAFTENSMUST00000075387JUNC00054', 'CGACAFTENST00000409762JUNC00056', 'CGACAFTENST00000295894JUNC00059', 'CGACAFTENST00000330622JUNC00068', 'CGACAFTENST00000456095JUNC00069', 'CGACAFTENST00000344384JUNC00075', 'CGACAFTENST00000357250JUNC00082', 'CGACAFTENST00000357250JUNC00083', 'CGACAFTENST00000357250JUNC00084', 'CGACAFTENST00000357250JUNC00085', 'CGACAFTENST00000320755JUNC00086', 'CGACAFTENST00000222597JUNC00087', 'CGACAFTENST00000295987JUNC00089', 'CGACAFTENST00000373836JUNC00091', 'CGACAFTENST00000373701JUNC00093', 'CGACAFTENST00000353341JUNC00094', 'CGACAFTENSMUST00000026416JUNC00100', 'CGACAFTENSMUST00000026416JUNC00101', 'CGACAFTENST00000315266JUNC00102', 'CGACAFTENSMUST00000047645JUNC00105', 'CGACAFTENST00000554684JUNC00106', 'CGACAFTENST00000356644JUNC00111', 'CGACAFTENST00000398376JUNC00112', 'CGACAFTENST00000398376JUNC00113', 'CGACAFTENSMUST00000160629JUNC00119', 'CGACAFTENST00000444053JUNC00120', 'CGACAFTENST00000630269JUNC00121', 'CGACAFTENST00000396840JUNC00131', 'CGACAFTENST00000451813JUNC00133', 'CGACAFTENST00000402711JUNC00134', 'CGACAFTENSMUST00000096766JUNC00135', 'CGACAFTENSMUST00000096766JUNC00136', 'CGACAFTENSMUST00000080934JUNC00138', 'CGACAFTENSMUST00000108191JUNC00139', 'CGACAFTENST00000379698JUNC00145', 'CGACAFTENST00000257468JUNC00149', 'CGACAFTENST00000263280JUNC00150', 'CGACAFTENST00000269137JUNC00151', 'CGACAFTENST00000355665JUNC00153', 'CGACAFTENST00000579690JUNC00154', 'CGACAFTENST00000444564JUNC00155', 'CGACAFTENST00000354600JUNC00159', 'CGACAFTENST00000299259JUNC00160', 'CGACAFTENST00000342556JUNC00161', 'CGACAFTENST00000345807JUNC00163', 'CGACAFTENST00000532317JUNC00168', 'CGACAFTENST00000393346JUNC00169', 'CGACAFTENSMUST00000207225JUNC00170', 'CGACAFTENSMUST00000114859JUNC00175', 'CGACAFTENST00000515064JUNC00177', 'CGACAFTENST00000515064JUNC00178', 'CGACAFTENST00000502640JUNC00180', 'CGACAFTENSMUST00000199256JUNC00181', 'CGACAFTENST00000313478JUNC00182', 'CGACAFTENST00000392020JUNC00183', 'CGACAFTENSMUST00000113190JUNC00185', 'CGACAFTENSMUST00000113190JUNC00186', 'CGACAFTENST00000373069JUNC00207', 'CGACAFTENST00000373069JUNC00208', 'CGACAFTENST00000373066JUNC00209', 'CGACAFTENST00000373064JUNC00210', 'CGACAFTENST00000343315JUNC00214', 'CGACAFTENST00000266732JUNC00215', 'CGACAFTENST00000393053JUNC00216', 'CGACAFTENSMUST00000099355JUNC00217', 'CGACAFTENST00000586534JUNC00227', 'CGACAFTENST00000391851JUNC00234', 'CGACAFTENST00000357167JUNC00240', 'CGACAFTENST00000218008JUNC00241', 'CGACAFTENST00000252485JUNC00242', 'CGACAFTENST00000544474JUNC00243', 'CGACAFTENST00000454580JUNC00245']
list_83_via_tr = list()

### TEST DES 83
nb_test_83 = 0
list_test_83 = list()
with open('/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/junction_validation_dog_2019-10-02/analysed_specific_junction_not_in_feelnc/specific_junction_83.bed', 'r') as test_file_83_i:
    lecture = test_file_83_i.readlines()
    for line in lecture:
        line = line.replace('\n', '').split('\t')
        tr = line[3].split('JUNC')[0]
        nb_test_83 += 1
        if not tr in list_test_83:
            list_test_83.append(tr)

print("Nb junction test:", nb_test_83)
print("Nb tr test:", len(list_test_83), '\n')
        
#for junction in list_83:
#    junction = junction.replace('\n', '').split('JUNC')[0]
#    if not junction in list_83_via_tr:
#        list_83_via_tr.append(junction)
##print("Nb tr ref:", len(list_83_via_tr))
##print()

### TEST DES 903
with open('/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/junction_validation_dog_2019-10-02/analysed_specific_junction_not_in_feelnc/specific_junction_cga_via_feelnc.bed', 'r') as ref_file_i:
    """
    ICI LES 903
    """
    lecture = ref_file_i.readlines()
    for line in lecture:
        line = line.replace('\n', '').split('\t')
        tr = line[3][:-9]
        
        if not (line[0], line[1], line[2], line[5], line[9], line[10], line[11]) in list_ref:
            list_ref.append((line[0], line[1], line[2], line[5], line[9], line[10], line[11]))
            nb_test += 1
        if tr not in list_ref_tr:
                list_ref_tr.append(tr)
            
        if tr not in dict_ref:
            dict_ref[tr] = []
        dict_ref[tr].append((line[0], line[1], line[2], line[5], line[8], line[9], line[10], line[11]))
        nb_junction_ref += 1
        
#    print("Nb tr with spe junctions in ref:", len(dic_ref))
    print("Nb junction in list ref:", len(list_ref))
    print("Nb junction in ref:", nb_junction_ref)
    print("Nb tr in ref:", len(list_ref_tr))
    print("Nb test:", nb_test)
    print()

nb_test = 0

### TEST DES 245
with open('/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/junction_validation_dog_2019-10-02/analysed_specific_junction_not_in_feelnc/specific_junction_cga2019.bed', 'r') as pred_file_i:
    """
    ICI LES 245
    """
    lecture = pred_file_i.readlines()
    for line in lecture:
        line = line.replace('\n', '').split('\t')
        tr = line[3][:-9]
        if not (line[0], line[1], line[2], line[5], line[9], line[10], line[11]) in list_pred:
            if (line[0], line[1], line[2], line[5], line[9], line[10], line[11]) in list_ref:
                list_pred.append((line[0], line[1], line[2], line[5], line[9], line[10], line[11]))
                nb_test += 1
        if tr not in list_pred_tr:
            list_pred_tr.append(tr)
        
#        if (line[0], line[1], line[2], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11]) in list_ref:
#            nb_test += 1 
        
        if tr not in dict_pred:
            dict_pred[tr] = []
        dict_pred[tr].append((line[0], line[1], line[2], line[5], line[9], line[10], line[11]))
        nb_junction_pred += 1

#    print("Nb tr with spe junctions in pred:", len(dict_pred))
    print("Nb junction in list pred:", len(list_pred))
    print("Nb junction in pred:", nb_junction_pred)
    print("Nb tr in pred:", len(list_pred_tr))
    print("Nb test:", nb_test)
    print()


nb_test = 0


print(len(list_ref_tr))
print(len(list_pred_tr))
with open("test.txt", 'w') as test_o:
    for junction in list_pred_tr:  
#        print(junction)
        if junction in list_83_via_tr:
            nb_test += 1
            test_o.write(str(junction)+'\n')
print(nb_test)


list_111 = ['CGACAFTENST00000257468', 'CGACAFTENSMUST00000047645', 'CGACAFTENST00000391851', 'CGACAFTENST00000263280', 'CGACAFTENST00000252485', 'CGACAFTENST00000454580', 'CGACAFTENSMUST00000026416', 'CGACAFTENST00000342556', 'CGACAFTENST00000357167', 'CGACAFTENST00000403231', 'CGACAFTENST00000630269', 'CGACAFTENST00000380036', 'CGACAFTENST00000297620', 'CGACAFTENSMUST00000114859', 'CGACAFTENST00000320755', 'CGACAFTENST00000357250', 'CGACAFTENSMUST00000080934', 'CGACAFTENSMUST00000108191', 'CGACAFTENST00000356644', 'CGACAFTENST00000526683', 'CGACAFTENST00000456095', 'CGACAFTENST00000527197', 'CGACAFTENST00000579690', 'CGACAFTENST00000355665', 'CGACAFTENST00000346177', 'CGACAFTENST00000343315', 'CGACAFTENSMUST00000099355', 'CGACAFTENST00000393053', 'CGACAFTENST00000266732', 'CGACAFTENSMUST00000096766', 'CGACAFTENST00000402711', 'CGACAFTENST00000444564', 'CGACAFTENST00000409744', 'CGACAFTENST00000410020', 'CGACAFTENST00000410041', 'CGACAFTENST00000394120', 'CGACAFTENST00000413539', 'CGACAFTENST00000258104', 'CGACAFTENST00000409762', 'CGACAFTENST00000360244', 'CGACAFTENST00000222597', 'CGACAFTENST00000448981', 'CGACAFTENST00000542836', 'CGACAFTENST00000373836', 'CGACAFTENSMUST00000105847', 'CGACAFTENST00000295894', 'CGACAFTENST00000396969', 'CGACAFTENST00000586534', 'CGACAFTENST00000393346', 'CGACAFTENSMUST00000207484', 'CGACAFTENSMUST00000207225', 'CGACAFTENST00000532317', 'CGACAFTENSMUST00000208730', 'CGACAFTENST00000379472', 'CGACAFTENST00000457054', 'CGACAFTENST00000451813', 'CGACAFTENST00000613795', 'CGACAFTENSMUST00000113190', 'CGACAFTENST00000392020', 'CGACAFTENST00000330622', 'CGACAFTENSMUST00000075387', 'CGACAFTENST00000404234', 'CGACAFTENSMUST00000100772', 'CGACAFTENST00000396840', 'CGACAFTENST00000380138', 'CGACAFTENST00000514493', 'CGACAFTENST00000615665', 'CGACAFTENST00000318445', 'CGACAFTENST00000488267', 'CGACAFTENST00000502640', 'CGACAFTENST00000515064', 'CGACAFTENSMUST00000199256', 'CGACAFTENST00000304523', 'CGACAFTENST00000299259', 'CGACAFTENST00000313478', 'CGACAFTENSMUST00000072750', 'CGACAFTENST00000353341', 'CGACAFTENST00000309703', 'CGACAFTENST00000354600', 'CGACAFTENST00000392456', 'CGACAFTENST00000444053', 'CGACAFTENST00000345807', 'CGACAFTENST00000398376', 'CGACAFTENST00000542447', 'CGACAFTENST00000537058', 'CGACAFTENSMUST00000114548', 'CGACAFTENST00000536727', 'CGACAFTENSMUST00000171041', 'CGACAFTENST00000377298', 'CGACAFTENSMUST00000160629', 'CGACAFTENST00000436336', 'CGACAFTENST00000370197', 'CGACAFTENST00000426398', 'CGACAFTENST00000260506', 'CGACAFTENSMUST00000174397', 'CGACAFTENST00000269137', 'CGACAFTENST00000315266', 'CGACAFTENST00000478722', 'CGACAFTENST00000554684', 'CGACAFTENST00000373064', 'CGACAFTENST00000432073', 'CGACAFTENST00000373066', 'CGACAFTENST00000373069', 'CGACAFTENST00000344384', 'CGACAFTENST00000379698', 'CGACAFTENSMUST00000115438', 'CGACAFTENST00000442742', 'CGACAFTENST00000295987', 'CGACAFTENST00000373701', 'CGACAFTENST00000218008', 'CGACAFTENST00000544474']
print()
print(len(list_test_83))
print(len(list_ref_tr))
print(len(list_pred_tr))

nb_test = 0
for test_83 in list_test_83:
    if test_83 in list_pred_tr:
        nb_test += 1
        
print("Test result: 74 in 208:", len(list_test_83), "in", len(list_pred_tr), ":", nb_test)

nb_test = 0
for test_83 in list_test_83:
    if test_83 in list_ref_tr:
        nb_test += 1
        
print("Test result: 74 in 166:", len(list_test_83), "in", len(list_ref_tr), ":", nb_test)
     
nb_test = 0
for test_83 in list_test_83:
    if test_83 in list_ref_tr and test_83 in list_pred_tr:
        nb_test += 1
        
print("Test result: 74 in 166 and 208:", len(list_test_83), "in", len(list_ref_tr), "and", len(list_pred_tr), ":", nb_test)





nb_test = 0
for test_111 in list_111:
    if test_111 in list_ref_tr and test_111 in list_pred_tr:
        nb_test += 1
        
print("Test result: 111 in 166 and 208:", len(list_111), "in", len(list_ref_tr), "and", len(list_pred_tr), ":", nb_test)

nb_test = 0
for test_111 in list_111:
    if test_111 in list_ref_tr:
        nb_test += 1
        
print("Test result: 111 in 166:", len(list_111), "in", len(list_ref_tr), ":", nb_test)

nb_test = 0
for test_111 in list_111:
    if test_111 in list_pred_tr:
        nb_test += 1
        
print("Test result: 111 in 208:", len(list_111), "in", len(list_pred_tr), ":", nb_test)

nb_test = 0
for tr_pred in list_pred_tr:
    if tr_pred in list_ref_tr and tr_pred not in list_test_83:
        nb_test += 1
        
print("Test result: 208 in 166:", nb_test)

nb_test = 0
for tr_111 in list_111:
    if tr_111 in list_test_83:
        nb_test +=1
        print(tr_111)
        
print("Nb 74 in 111:", nb_test, len(list_test_83))

with open('/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/junction_validation_dog_2019-10-02/analysed_specific_junction_not_in_feelnc/specific_junction_83.bed', 'r') as test_file_83_i:
    lecture = test_file_83_i.readlines()
    for line in lecture:
        line = line.replace('\n', '').split('\t')
        tr = line[3].split('JUNC')[0]
        if tr == "CGACAFTENST00000257468":
            print(line)

list_pred_tr
list_ref_tr

#test_nb_pred = 0
#test_nb_ref = 0
#
#for tr_pred in dict_pred:
#    test_nb_pred += len(dict_pred[tr_pred])
#
#for tr_ref in dict_ref:
#    test_nb_ref += len(dict_ref[tr_ref])
#    
#print(">>> PREDICTION <<<")
#print("nb_junctions:", test_nb_pred)
#print("nb_transcripts:", len(dict_pred))
#print(">>> REFERENCE <<<")
#print("nb_junctions:", test_nb_ref)
#print("nb_transcripts:", len(dict_ref))
#
#dict_conserved_prediction_specificity = dict()
#
#test_nb_in_ref = 0
#test_nb_not_in_ref = 0
#
#junc = 0
#
#with open("specific_junction_cga_not_in_feelnc.bed", 'w') as output_s:
#    
#    for tr_pred in dict_pred:
#        for junction_pred in dict_pred[tr_pred]:
#            junc += 1
#            try:
#                if junction_pred in dict_ref[tr_pred]:
#                    test_nb_in_ref += 1
#                else:
#                    test_nb_not_in_ref += 1
#                    if junc < 10:
#                        output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC0000"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                    elif junc < 100 :
#                        output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC000"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                    elif junc < 1000:
#                        output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC00"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                    elif junc < 10000:
#                        output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC0"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                    else:
#                        output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#
#            except KeyError:
#                test_nb_not_in_ref +=1
#                if junc < 10:
#                    output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC0000"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                elif junc < 100 :
#                    output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC000"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                elif junc < 1000:
#                    output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC00"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                elif junc < 10000:
#                    output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC0"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#                else:
#                    output_s.write(junction_pred[0] +'\t'+ junction_pred[1] +'\t'+ junction_pred[2] +'\t'+ tr_pred+"JUNC"+str(junc) +'\t'+ junction_pred[3] +'\t'+ junction_pred[4] +'\t'+ junction_pred[5] +'\t'+ junction_pred[6] +'\t'+ junction_pred[7] +'\t'+ junction_pred[8] +'\t'+ junction_pred[9] +'\t'+ junction_pred[10]+'\n')
#
#        
#print("nb junctions in ref:", test_nb_in_ref)
#print("nb junction not in ref", test_nb_not_in_ref)
#
#with open("list_ref.txt", 'w') as output_r:
#    for tr in dict_ref.values():
#        for tr_i in tr:
#            output_r.write(str(tr_i)+'\n')
#            
#with open("list_pred.txt", 'w') as output_p:
#    for tr in dict_pred.values():
#        for tr_i in tr:
#            output_p.write(str(tr_i)+'\n')