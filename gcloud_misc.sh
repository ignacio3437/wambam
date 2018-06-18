#Things about gcloud


#To undo agalma messing up the python path
export BOTO_CONFIG=/dev/null


#To copy files from your bucket to the cloud
gsutil cp gcs://bucket-name/you-file .


#SSH behind a proxy
ssh USER@FINAL_DEST -o "ProxyCommand=nc -X connect -x PROXYHOST:PROXYPORT %h %p"

#ssh to glcloud
ssh -i ~/.ssh/gc_rsa ignacio3437@35.198.226.23 -o "ProxyCommand=nc -X connect -x 10.150.2.32:80 %h %p"

#to generate key
ssh-keygen -t rsa -f ~/.ssh/gc_rsa -C ignacio3437 
#Pas for this is a*s

