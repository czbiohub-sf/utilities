MULTIOMICS_ECR=423543210473.dkr.ecr.us-west-2.amazonaws.com/multiomics

multiomics/cellranger-arc-2.0.0.tar.gz:
	$(error You need to obtain a copy of $@ and place it in multiomics/$@)

multiomics/cellranger-6.0.2.tar.gz:
	$(error You need to obtain a copy of $@ and place it in multiomics/$@)

build-multiomics-docker: multiomics/cellranger-arc-2.0.0.tar.gz multiomics/cellranger-6.0.2.tar.gz
	cd multiomics && \
	docker build . --no-cache -t multiomics:latest

push-multiomics-docker:
	docker tag multiomics:latest $(MULTIOMICS_ECR):latest
	docker push $(MULTIOMICS_ECR)

release-multiomics-docker: build-multiomics-docker push-multiomics-docker