# platelet-wcm webapp image
# Extends the runtime image with project source and the Dash web UI entrypoint.
#
# Build locally from the platelet-wcm/ project root:
#
#     docker build -f docker/webapp/Dockerfile \
#         --build-arg from=platelet-wcm-runtime \
#         -t platelet-wcm-webapp .
#
# Run locally:
#
#     docker run --rm -p 8050:8050 -v $(pwd)/out:/platelet-wcm/out \
#         -e WEBAPP_PORT=8050 platelet-wcm-webapp
#
#     Then open http://localhost:8050/

ARG from=platelet-wcm-runtime:latest
FROM ${from}

ARG git_hash=""
ARG git_branch=""
ENV IMAGE_GIT_HASH="$git_hash" \
    IMAGE_GIT_BRANCH="$git_branch"

LABEL application="Platelet Whole-Cell Model" \
    repository="https://github.com/stevehaigh/platelet-wcm"

COPY . /platelet-wcm
WORKDIR /platelet-wcm

ENV PYTHONPATH=/platelet-wcm
ENV WCECOLI_WEBAPP_MODE=container
ENV WEBAPP_PORT=80
EXPOSE 80

CMD ["sh", "-c", "python runscripts/manual/webapp.py --port ${WEBAPP_PORT}"]
