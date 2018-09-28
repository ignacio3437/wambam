#Things about gcloud


#To undo agalma messing up the python path
export BOTO_CONFIG=/dev/null


#To copy files to bucket
gsutil cp file.txt gs://bucket-name

#To copy files from your bucket to the cloud
gsutil cp gcs://bucket-name/you-file .
gsutil cp gs://pteraeolidia1/Coral_public_txtm.zip .



#Fireup Ubuntu Instances

#to generate key
ssh-keygen -t rsa -f ~/.ssh/gc_rsa_wammsu -C wammsu 



#SSH behind a proxy
# ssh USER@FINAL_DEST -o "ProxyCommand=nc -X connect -x PROXYHOST:PROXYPORT %h %p"

#ssh to glcloud
ssh -i ~/.ssh/gc_rsa_wammsu wammsu@35.234.0.89 -o "ProxyCommand=nc -X connect -x 10.150.2.32:80 %h %p"

#known host error
rm ~/.ssh/known_hosts

