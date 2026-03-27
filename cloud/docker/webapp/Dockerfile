# Container image #3: wcm-webapp.
# Extends wcm-code with the Dash web UI as the entrypoint.
#
# To build locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/webapp/Dockerfile \
#         --build-arg from=${USER}-wcm-code \
#         -t ${USER}-wcm-webapp .
#
# To run locally:
#
#     > docker run --rm -p 8050:8050 -v $(pwd)/out:/wcEcoli/out ${USER}-wcm-webapp
#
# Then open http://localhost:8050/
#
# In the cloud (GitHub Actions deploys to Azure ACI on port 80):
#
#     http://wcecoli-webapp.uksouth.azurecontainer.io/

ARG from=wcm-code:latest
FROM ${from}

# Copy webapp source files (not present in the base wcm-code image)
COPY wholecell/webapp/ /wcEcoli/wholecell/webapp/
COPY runscripts/manual/webapp.py /wcEcoli/runscripts/manual/webapp.py

# The webapp uses local subprocess for simulations (no docker-in-docker).
ENV WCECOLI_WEBAPP_MODE=container

# Default to port 80 for cloud; override with --port for local dev.
ENV WEBAPP_PORT=80
EXPOSE 80

CMD ["sh", "-c", "python runscripts/manual/webapp.py --port ${WEBAPP_PORT}"]
