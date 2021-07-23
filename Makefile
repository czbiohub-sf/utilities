MULTIOMICS_ECR=423543210473.dkr.ecr.us-west-2.amazonaws.com/multiomics

multiomics/cellranger-arc-2.0.0.tar.gz:
	$(error You need to obtain a copy of $@ and place it here)

build-multiomics-docker: multiomics/cellranger-arc-2.0.0.tar.gz
	cd multiomics && \
	docker build . -t multiomics:latest

push-multiomics-docker: build-multiomics-docker
	docker tag multiomics:latest $(MULTIOMICS_ECR):latest
	docker push $(MULTIOMICS_ECR)
