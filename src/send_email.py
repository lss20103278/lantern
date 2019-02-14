#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: email.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Mon 03 Dec 2018 03:25:19 PM CST
#########################################################################

import sys
import os

kwargs_raw = sys.argv[1:]
kwargs={'-send':''}
for i in range(len(kwargs_raw)):
    for j in kwargs.keys():
        if(kwargs_raw[i]==j):
            kwargs.update({j:kwargs_raw[i+1]})
choose_send = kwargs.get('-send')

dirname = '/'.join(os.getcwd().split('/')[-2:])
for f in os.listdir('.'):
    if 'eml' in f:
        sn = f.split(']')[1].split('.eml')[0]
        #sn = f[1:-18]
    if 'zip' in f:
        filename = f

#ofile = open('email.txt','w')
#ofile.write('主题：\n'+sn+'结果\n')
#ofile.write('内容：\n以下为儿童遗传病分析结果：\n批次：'+sn+'结果\n\n位置\n/DATA/sslyu/Project/Genetics_children_hospital/2018/'+dirname+'\n\n祝好，\n吕珊珊\n\n')
#ofile.close()

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

email_user='sslv@basepair.cn'
email_send1='lss@sibs.ac.cn'
email_send2=['yyzhou@basepair.cn', 'jjia@basepair.cn', 'gli@basepair.cn', 'yjsun@basepair.cn'] # in smtplib module, when multi recipients, need a list 
#email_send2=['yyzhou@basepair.cn', 'jjia@basepair.cn', 'gli@basepair.cn', 'yjsun@basepair.cn', 'jylou@basepair.cn'] # in smtplib module, when multi recipients, need a list 

subject = sn+'结果'

msg = MIMEMultipart()
msg['From'] = email_user
#choose_send = raw_input('please choose 1(lss@sibs.ac.cn) or 2(yyzhou@basepair.cn,jjia@basepair.cn,gli@basepair.cn,yjsun@basepair.cn):  ')
if choose_send == '1':
    msg['To'] = email_send1
else:
    msg['To'] = ", ".join(email_send2) # in email module, need a string
msg['Subject'] = subject

if os.path.exists('email.txt'):
    body = open('email.txt').readlines()
    body = ''.join(body)
else:
    print """
    Run /DATA/sslyu/trio_BWA-GATK_3.0/src/create_email_content.py first
    """
    sys.exit()

#msg.attach(MIMEText(body, 'plain')) # can't display Chineses character 
#msg.attach(MIMEText(body, 'html', 'utf-8')) # format is problematic
msg.attach(MIMEText(body, 'plain', 'utf-8'))
attachment = open(filename, 'rb')

part = MIMEBase('application', 'octet-stream')
part.set_payload((attachment).read())
encoders.encode_base64(part)
part.add_header('Content-Disposition', "attachment; filename="+filename)

msg.attach(part)
text = msg.as_string()

mail = smtplib.SMTP('smtp.basepair.cn', 25)
mail.ehlo()
mail.starttls()
mail.login(email_user, 'azby06078961!')

if choose_send == '1':
    mail.sendmail(email_user, email_send1, text)
else:
    mail.sendmail(email_user, email_send2, text)

mail.quit()
