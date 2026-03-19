# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import tempfile
import yaml

from lsst.daf.butler import Butler, Config
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.pipe.base import Pipeline

with open("metricDescriptions.yaml") as f:
    metricDescriptions = yaml.safe_load(f)

class NoAliasDumper(yaml.Dumper):
    def ignore_aliases(self, data):
        return True


class QuotedStr(str):
    pass


def quoted_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style="'")


NoAliasDumper.add_representer(QuotedStr, quoted_representer)


DRP_PIPE_DIR = os.environ['DRP_PIPE_DIR']
DRP_FILE = os.path.join(DRP_PIPE_DIR, 'pipelines/HSC/DRP-RC2.yaml')
TEST_DIR = os.path.abspath(os.path.dirname(__file__))

root = makeTestTempDir(TEST_DIR)
tmpdir = tempfile.mkdtemp(dir=root)

config = Config()
config["registry", "db"] = f"sqlite:///{tmpdir}/gen3.sqlite3"
config = Butler.makeRepo(root, config)
butler = Butler.from_config(config, writeable=True)
pipeline = Pipeline.from_uri(DRP_FILE)
pipeline_graph = pipeline.to_graph(registry=butler.registry,
                                   visualization_only=True)

# Generate a lookup table for bundles that are put into tables:
tableNameLookup = {}
for taskName in pipeline_graph.tasks:
    if 'makeMetricTable' in pipeline_graph.tasks[taskName].task_class_name:
        task = pipeline_graph.tasks[taskName]
        inputBundle = task.inputs['data'].dataset_type_name
        outputTable = task.outputs['metricTable'].dataset_type_name
        tableNameLookup[inputBundle] = outputTable

metricData = {
    'lowThreshold': None,
    'highThreshold': None,
    'metricTypes': None,
    'description': None,
}

metricInfo = {'tasks': {}}
for taskName in pipeline_graph.tasks:
    task = pipeline_graph.tasks[taskName]
    if 'metrics' in task.outputs.keys():
        metricBundleName = task.outputs['metrics'].dataset_type_name
        taskConfigDict = task.config.toDict()
        atools = taskConfigDict['atools']
        tableName = tableNameLookup.get(metricBundleName)
        for atool in atools:
            if 'metric' in taskConfigDict['atools'][atool]['produce']:
                metricProduce = taskConfigDict['atools'][atool]['produce']['metric']
                if len(metrics := metricProduce['units'].keys()) > 0:
                    print(metricBundleName)
                    metricInfo['tasks'].setdefault(taskName, {
                        'bundleName': metricBundleName,
                        'tableName': tableName,
                        'atools': {}
                    })
                    metricInfo['tasks'][taskName]['atools'].setdefault(atool, {})
                    for metric in metrics:
                        metricName = metricProduce['newNames'].get(metric, metric)
                        if metricName[0] == '{' or ' ' in metricName:
                            metricName = QuotedStr(metricName)
                        descriptionData = (metricDescriptions
                                           .get('tasks', {})
                                           .get(taskName, {})
                                           .get('atools', {})
                                           .get(atool, {})
                                           .get(metricName, {}))
                        entry = metricData.copy()
                        entry['metricTypes'] = descriptionData.get('metricTypes')
                        entry['description'] = descriptionData.get('description')
                        metricInfo['tasks'][taskName]['atools'][atool] |= {metricName: entry}

with open("output.yaml", "w") as f:
    yaml.dump(metricInfo, f,
              Dumper=NoAliasDumper,
              sort_keys=False,
              default_flow_style=False,
              indent=2)

removeTestTempDir(root)