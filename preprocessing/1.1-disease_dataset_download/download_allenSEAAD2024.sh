# path to git repository: https://github.com/AllenInstitute/SEA-AD_2024
# path to AWS bucket: https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#PFC/RNAseq/

# install aws cli
#curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
#unzip awscliv2.zip
#sudo ./aws/install

# Update pyopenssl
#python3 -m venv /home/gugene/venv/
#/home/gugene/venv/bin/pip install pyopenssl --upgrade

cd ${DATA_DIR:-${DATA_DIR:-/mnt/bdata/gugene}}/datasets/SN_RNAseq/sea_ad_2024
/home/gugene/aws/dist/aws s3 ls --no-sign-request s3://sea-ad-single-cell-profiling/PFC/RNAseq/
/home/gugene/aws/dist/aws s3 cp s3://sea-ad-single-cell-profiling/PFC/RNAseq/ RNAseq --recursive --no-sign-request
