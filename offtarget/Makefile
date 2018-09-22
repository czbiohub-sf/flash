fmt:
	go fmt .
vet: fmt
	go vet .

test: fmt vet
	go test . -cover

binary:
	GOOS=linux go build .

image: binary
	docker build -t ryansking/crispr_offtarget:latest .

.PHONY: fmt vet test image binary
