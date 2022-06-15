/*
 * Copyright (c) 2019 TAOS Data, Inc. <jhtao@taosdata.com>
 *
 * This program is free software: you can use, redistribute, and/or modify
 * it under the terms of the GNU Affero General Public License, version 3
 * or later ("AGPL"), as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "functionMgt.h"
#include "os.h"
#include "parAst.h"
#include "parInt.h"
#include "parToken.h"
#include "systable.h"
#include "tglobal.h"

typedef void* (*FMalloc)(size_t);
typedef void (*FFree)(void*);

extern void* ParseAlloc(FMalloc);
extern void  Parse(void*, int, SToken, void*);
extern void  ParseFree(void*, FFree);
extern void  ParseTrace(FILE*, char*);

int32_t parse(SParseContext* pParseCxt, SQuery** pQuery) {
  SAstCreateContext cxt;
  initAstCreateContext(pParseCxt, &cxt);
  void*   pParser = ParseAlloc((FMalloc)taosMemoryMalloc);
  int32_t i = 0;
  while (1) {
    SToken t0 = {0};
    if (cxt.pQueryCxt->pSql[i] == 0) {
      Parse(pParser, 0, t0, &cxt);
      goto abort_parse;
    }
    t0.n = tGetToken((char*)&cxt.pQueryCxt->pSql[i], &t0.type);
    t0.z = (char*)(cxt.pQueryCxt->pSql + i);
    i += t0.n;

    switch (t0.type) {
      case TK_NK_SPACE:
      case TK_NK_COMMENT: {
        break;
      }
      case TK_NK_SEMI: {
        Parse(pParser, 0, t0, &cxt);
        goto abort_parse;
      }
      case TK_NK_ILLEGAL: {
        snprintf(cxt.pQueryCxt->pMsg, cxt.pQueryCxt->msgLen, "unrecognized token: \"%s\"", t0.z);
        cxt.errCode = TSDB_CODE_PAR_SYNTAX_ERROR;
        goto abort_parse;
      }
      case TK_NK_HEX:
      case TK_NK_OCT:
      case TK_NK_BIN: {
        snprintf(cxt.pQueryCxt->pMsg, cxt.pQueryCxt->msgLen, "unsupported token: \"%s\"", t0.z);
        cxt.errCode = TSDB_CODE_PAR_SYNTAX_ERROR;
        goto abort_parse;
      }
      default:
        // ParseTrace(stdout, "");
        Parse(pParser, t0.type, t0, &cxt);
        if (TSDB_CODE_SUCCESS != cxt.errCode) {
          goto abort_parse;
        }
    }
  }

abort_parse:
  ParseFree(pParser, (FFree)taosMemoryFree);
  if (TSDB_CODE_SUCCESS == cxt.errCode) {
    *pQuery = (SQuery*)nodesMakeNode(QUERY_NODE_QUERY);
    if (NULL == *pQuery) {
      return TSDB_CODE_OUT_OF_MEMORY;
    }
    (*pQuery)->pRoot = cxt.pRootNode;
    (*pQuery)->placeholderNum = cxt.placeholderNo;
    TSWAP((*pQuery)->pPlaceholderValues, cxt.pPlaceholderValues);
  }
  taosArrayDestroy(cxt.pPlaceholderValues);
  return cxt.errCode;
}

typedef struct SCollectMetaKeyCxt {
  SParseContext*   pParseCxt;
  SParseMetaCache* pMetaCache;
  SNode*           pStmt;
} SCollectMetaKeyCxt;

typedef struct SCollectMetaKeyFromExprCxt {
  SCollectMetaKeyCxt* pComCxt;
  int32_t             errCode;
} SCollectMetaKeyFromExprCxt;

static int32_t collectMetaKeyFromQuery(SCollectMetaKeyCxt* pCxt, SNode* pStmt);

static EDealRes collectMetaKeyFromFunction(SCollectMetaKeyFromExprCxt* pCxt, SFunctionNode* pFunc) {
  if (fmIsBuiltinFunc(pFunc->functionName)) {
    return DEAL_RES_CONTINUE;
  }
  pCxt->errCode = reserveUdfInCache(pFunc->functionName, pCxt->pComCxt->pMetaCache);
  return TSDB_CODE_SUCCESS == pCxt->errCode ? DEAL_RES_CONTINUE : DEAL_RES_ERROR;
}

static bool needGetTableIndex(SNode* pStmt) {
  if (QUERY_SMA_OPTIMIZE_ENABLE == tsQuerySmaOptimize && QUERY_NODE_SELECT_STMT == nodeType(pStmt)) {
    SSelectStmt* pSelect = (SSelectStmt*)pStmt;
    return (NULL != pSelect->pWindow && QUERY_NODE_INTERVAL_WINDOW == nodeType(pSelect->pWindow));
  }
  return false;
}

static int32_t collectMetaKeyFromRealTableImpl(SCollectMetaKeyCxt* pCxt, SRealTableNode* pRealTable,
                                               AUTH_TYPE authType) {
  int32_t code = reserveTableMetaInCache(pCxt->pParseCxt->acctId, pRealTable->table.dbName, pRealTable->table.tableName,
                                         pCxt->pMetaCache);
  if (TSDB_CODE_SUCCESS == code) {
    code = reserveTableVgroupInCache(pCxt->pParseCxt->acctId, pRealTable->table.dbName, pRealTable->table.tableName,
                                     pCxt->pMetaCache);
  }
  if (TSDB_CODE_SUCCESS == code) {
    code = reserveUserAuthInCache(pCxt->pParseCxt->acctId, pCxt->pParseCxt->pUser, pRealTable->table.dbName, authType,
                                  pCxt->pMetaCache);
  }
  if (TSDB_CODE_SUCCESS == code) {
    code = reserveDbVgInfoInCache(pCxt->pParseCxt->acctId, pRealTable->table.dbName, pCxt->pMetaCache);
  }
  if (TSDB_CODE_SUCCESS == code && needGetTableIndex(pCxt->pStmt)) {
    code = reserveTableIndexInCache(pCxt->pParseCxt->acctId, pRealTable->table.dbName, pRealTable->table.tableName,
                                    pCxt->pMetaCache);
  }
  return code;
}

static EDealRes collectMetaKeyFromRealTable(SCollectMetaKeyFromExprCxt* pCxt, SRealTableNode* pRealTable) {
  pCxt->errCode = collectMetaKeyFromRealTableImpl(pCxt->pComCxt, pRealTable, AUTH_TYPE_READ);
  return TSDB_CODE_SUCCESS == pCxt->errCode ? DEAL_RES_CONTINUE : DEAL_RES_ERROR;
}

static EDealRes collectMetaKeyFromTempTable(SCollectMetaKeyFromExprCxt* pCxt, STempTableNode* pTempTable) {
  pCxt->errCode = collectMetaKeyFromQuery(pCxt->pComCxt, pTempTable->pSubquery);
  return TSDB_CODE_SUCCESS == pCxt->errCode ? DEAL_RES_CONTINUE : DEAL_RES_ERROR;
}

static EDealRes collectMetaKeyFromExprImpl(SNode* pNode, void* pContext) {
  SCollectMetaKeyFromExprCxt* pCxt = pContext;
  switch (nodeType(pNode)) {
    case QUERY_NODE_FUNCTION:
      return collectMetaKeyFromFunction(pCxt, (SFunctionNode*)pNode);
    case QUERY_NODE_REAL_TABLE:
      return collectMetaKeyFromRealTable(pCxt, (SRealTableNode*)pNode);
    case QUERY_NODE_TEMP_TABLE:
      return collectMetaKeyFromTempTable(pCxt, (STempTableNode*)pNode);
    default:
      break;
  }
  return DEAL_RES_CONTINUE;
}

static int32_t collectMetaKeyFromExprs(SCollectMetaKeyCxt* pCxt, SNodeList* pList) {
  SCollectMetaKeyFromExprCxt cxt = {.pComCxt = pCxt, .errCode = TSDB_CODE_SUCCESS};
  nodesWalkExprs(pList, collectMetaKeyFromExprImpl, &cxt);
  return cxt.errCode;
}

static int32_t collectMetaKeyFromSetOperator(SCollectMetaKeyCxt* pCxt, SSetOperator* pStmt) {
  int32_t code = collectMetaKeyFromQuery(pCxt, pStmt->pLeft);
  if (TSDB_CODE_SUCCESS == code) {
    code = collectMetaKeyFromQuery(pCxt, pStmt->pRight);
  }
  if (TSDB_CODE_SUCCESS == code) {
    code = collectMetaKeyFromExprs(pCxt, pStmt->pOrderByList);
  }
  return code;
}

static int32_t collectMetaKeyFromSelect(SCollectMetaKeyCxt* pCxt, SSelectStmt* pStmt) {
  SCollectMetaKeyFromExprCxt cxt = {.pComCxt = pCxt, .errCode = TSDB_CODE_SUCCESS};
  nodesWalkSelectStmt(pStmt, SQL_CLAUSE_FROM, collectMetaKeyFromExprImpl, &cxt);
  return cxt.errCode;
}

static int32_t collectMetaKeyFromAlterDatabase(SCollectMetaKeyCxt* pCxt, SAlterDatabaseStmt* pStmt) {
  return reserveDbCfgInCache(pCxt->pParseCxt->acctId, pStmt->dbName, pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromCreateTable(SCollectMetaKeyCxt* pCxt, SCreateTableStmt* pStmt) {
  int32_t code = reserveDbCfgInCache(pCxt->pParseCxt->acctId, pStmt->dbName, pCxt->pMetaCache);
  if (TSDB_CODE_SUCCESS == code && NULL == pStmt->pTags) {
    code = reserveTableVgroupInCache(pCxt->pParseCxt->acctId, pStmt->dbName, pStmt->tableName, pCxt->pMetaCache);
  }
  return code;
}

static int32_t collectMetaKeyFromCreateMultiTable(SCollectMetaKeyCxt* pCxt, SCreateMultiTableStmt* pStmt) {
  int32_t code = TSDB_CODE_SUCCESS;
  SNode*  pNode = NULL;
  FOREACH(pNode, pStmt->pSubTables) {
    SCreateSubTableClause* pClause = (SCreateSubTableClause*)pNode;
    code = reserveDbCfgInCache(pCxt->pParseCxt->acctId, pClause->dbName, pCxt->pMetaCache);
    if (TSDB_CODE_SUCCESS == code) {
      code =
          reserveTableMetaInCache(pCxt->pParseCxt->acctId, pClause->useDbName, pClause->useTableName, pCxt->pMetaCache);
    }
    if (TSDB_CODE_SUCCESS == code) {
      code = reserveTableVgroupInCache(pCxt->pParseCxt->acctId, pClause->dbName, pClause->tableName, pCxt->pMetaCache);
    }
    if (TSDB_CODE_SUCCESS != code) {
      break;
    }
  }
  return code;
}

static int32_t collectMetaKeyFromDropTable(SCollectMetaKeyCxt* pCxt, SDropTableStmt* pStmt) {
  int32_t code = TSDB_CODE_SUCCESS;
  SNode*  pNode = NULL;
  FOREACH(pNode, pStmt->pTables) {
    SDropTableClause* pClause = (SDropTableClause*)pNode;
    code = reserveTableMetaInCache(pCxt->pParseCxt->acctId, pClause->dbName, pClause->tableName, pCxt->pMetaCache);
    if (TSDB_CODE_SUCCESS == code) {
      code = reserveTableVgroupInCache(pCxt->pParseCxt->acctId, pClause->dbName, pClause->tableName, pCxt->pMetaCache);
    }
    if (TSDB_CODE_SUCCESS != code) {
      break;
    }
  }
  return code;
}

static int32_t collectMetaKeyFromAlterTable(SCollectMetaKeyCxt* pCxt, SAlterTableStmt* pStmt) {
  int32_t code = reserveTableMetaInCache(pCxt->pParseCxt->acctId, pStmt->dbName, pStmt->tableName, pCxt->pMetaCache);
  if (TSDB_CODE_SUCCESS == code) {
    code = reserveTableVgroupInCache(pCxt->pParseCxt->acctId, pStmt->dbName, pStmt->tableName, pCxt->pMetaCache);
  }
  return code;
}

static int32_t collectMetaKeyFromUseDatabase(SCollectMetaKeyCxt* pCxt, SUseDatabaseStmt* pStmt) {
  return reserveDbVgVersionInCache(pCxt->pParseCxt->acctId, pStmt->dbName, pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromCreateIndex(SCollectMetaKeyCxt* pCxt, SCreateIndexStmt* pStmt) {
  int32_t code = TSDB_CODE_SUCCESS;
  if (INDEX_TYPE_SMA == pStmt->indexType) {
    code = reserveTableMetaInCache(pCxt->pParseCxt->acctId, pCxt->pParseCxt->db, pStmt->tableName, pCxt->pMetaCache);
    if (TSDB_CODE_SUCCESS == code) {
      code =
          reserveTableVgroupInCache(pCxt->pParseCxt->acctId, pCxt->pParseCxt->db, pStmt->tableName, pCxt->pMetaCache);
    }
    if (TSDB_CODE_SUCCESS == code) {
      code = reserveDbVgInfoInCache(pCxt->pParseCxt->acctId, pCxt->pParseCxt->db, pCxt->pMetaCache);
    }
  }
  return code;
}

static int32_t collectMetaKeyFromCreateTopic(SCollectMetaKeyCxt* pCxt, SCreateTopicStmt* pStmt) {
  if (NULL != pStmt->pQuery) {
    return collectMetaKeyFromQuery(pCxt, pStmt->pQuery);
  }
  return TSDB_CODE_SUCCESS;
}

static int32_t collectMetaKeyFromExplain(SCollectMetaKeyCxt* pCxt, SExplainStmt* pStmt) {
  return collectMetaKeyFromQuery(pCxt, pStmt->pQuery);
}

static int32_t collectMetaKeyFromDescribe(SCollectMetaKeyCxt* pCxt, SDescribeStmt* pStmt) {
  SName name = {.type = TSDB_TABLE_NAME_T, .acctId = pCxt->pParseCxt->acctId};
  strcpy(name.dbname, pStmt->dbName);
  strcpy(name.tname, pStmt->tableName);
  int32_t code = catalogRemoveTableMeta(pCxt->pParseCxt->pCatalog, &name);
  if (TSDB_CODE_SUCCESS == code) {
    code = reserveTableMetaInCache(pCxt->pParseCxt->acctId, pStmt->dbName, pStmt->tableName, pCxt->pMetaCache);
  }
  return code;
}

static int32_t collectMetaKeyFromCreateStream(SCollectMetaKeyCxt* pCxt, SCreateStreamStmt* pStmt) {
  return collectMetaKeyFromQuery(pCxt, pStmt->pQuery);
}

static int32_t collectMetaKeyFromShowDnodes(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_DNODES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowMnodes(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_MNODES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowModules(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_MODULES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowQnodes(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_QNODES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowSnodes(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_SNODES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowBnodes(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_BNODES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowDatabases(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_USER_DATABASES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowFunctions(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_USER_FUNCTIONS,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowIndexes(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_USER_INDEXES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowStables(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_USER_STABLES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowStreams(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_PERFORMANCE_SCHEMA_DB, TSDB_PERFS_TABLE_STREAMS,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowTables(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  int32_t code = reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB,
                                         TSDB_INS_TABLE_USER_TABLES, pCxt->pMetaCache);
  if (TSDB_CODE_SUCCESS == code) {
    if (NULL != pStmt->pDbName) {
      code = reserveDbVgInfoInCache(pCxt->pParseCxt->acctId, ((SValueNode*)pStmt->pDbName)->literal, pCxt->pMetaCache);
    } else {
      code = reserveDbVgInfoInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, pCxt->pMetaCache);
    }
  }
  return code;
}

static int32_t collectMetaKeyFromShowUsers(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_USER_USERS,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowLicence(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_LICENCES,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowVgroups(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_INFORMATION_SCHEMA_DB, TSDB_INS_TABLE_VGROUPS,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowTopics(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_PERFORMANCE_SCHEMA_DB, TSDB_PERFS_TABLE_TOPICS,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromShowTransactions(SCollectMetaKeyCxt* pCxt, SShowStmt* pStmt) {
  return reserveTableMetaInCache(pCxt->pParseCxt->acctId, TSDB_PERFORMANCE_SCHEMA_DB, TSDB_PERFS_TABLE_TRANS,
                                 pCxt->pMetaCache);
}

static int32_t collectMetaKeyFromDelete(SCollectMetaKeyCxt* pCxt, SDeleteStmt* pStmt) {
  return collectMetaKeyFromRealTableImpl(pCxt, (SRealTableNode*)pStmt->pFromTable, AUTH_TYPE_WRITE);
}

static int32_t collectMetaKeyFromQuery(SCollectMetaKeyCxt* pCxt, SNode* pStmt) {
  pCxt->pStmt = pStmt;
  switch (nodeType(pStmt)) {
    case QUERY_NODE_SET_OPERATOR:
      return collectMetaKeyFromSetOperator(pCxt, (SSetOperator*)pStmt);
    case QUERY_NODE_SELECT_STMT:
      return collectMetaKeyFromSelect(pCxt, (SSelectStmt*)pStmt);
    case QUERY_NODE_ALTER_DATABASE_STMT:
      return collectMetaKeyFromAlterDatabase(pCxt, (SAlterDatabaseStmt*)pStmt);
    case QUERY_NODE_CREATE_TABLE_STMT:
      return collectMetaKeyFromCreateTable(pCxt, (SCreateTableStmt*)pStmt);
    case QUERY_NODE_CREATE_MULTI_TABLE_STMT:
      return collectMetaKeyFromCreateMultiTable(pCxt, (SCreateMultiTableStmt*)pStmt);
    case QUERY_NODE_DROP_TABLE_STMT:
      return collectMetaKeyFromDropTable(pCxt, (SDropTableStmt*)pStmt);
    case QUERY_NODE_ALTER_TABLE_STMT:
      return collectMetaKeyFromAlterTable(pCxt, (SAlterTableStmt*)pStmt);
    case QUERY_NODE_USE_DATABASE_STMT:
      return collectMetaKeyFromUseDatabase(pCxt, (SUseDatabaseStmt*)pStmt);
    case QUERY_NODE_CREATE_INDEX_STMT:
      return collectMetaKeyFromCreateIndex(pCxt, (SCreateIndexStmt*)pStmt);
    case QUERY_NODE_CREATE_TOPIC_STMT:
      return collectMetaKeyFromCreateTopic(pCxt, (SCreateTopicStmt*)pStmt);
    case QUERY_NODE_EXPLAIN_STMT:
      return collectMetaKeyFromExplain(pCxt, (SExplainStmt*)pStmt);
    case QUERY_NODE_DESCRIBE_STMT:
      return collectMetaKeyFromDescribe(pCxt, (SDescribeStmt*)pStmt);
    case QUERY_NODE_CREATE_STREAM_STMT:
      return collectMetaKeyFromCreateStream(pCxt, (SCreateStreamStmt*)pStmt);
    case QUERY_NODE_SHOW_DNODES_STMT:
      return collectMetaKeyFromShowDnodes(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_MNODES_STMT:
      return collectMetaKeyFromShowMnodes(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_MODULES_STMT:
      return collectMetaKeyFromShowModules(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_QNODES_STMT:
      return collectMetaKeyFromShowQnodes(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_SNODES_STMT:
      return collectMetaKeyFromShowSnodes(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_BNODES_STMT:
      return collectMetaKeyFromShowBnodes(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_DATABASES_STMT:
      return collectMetaKeyFromShowDatabases(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_FUNCTIONS_STMT:
      return collectMetaKeyFromShowFunctions(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_INDEXES_STMT:
      return collectMetaKeyFromShowIndexes(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_STABLES_STMT:
      return collectMetaKeyFromShowStables(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_STREAMS_STMT:
      return collectMetaKeyFromShowStreams(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_TABLES_STMT:
      return collectMetaKeyFromShowTables(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_USERS_STMT:
      return collectMetaKeyFromShowUsers(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_LICENCE_STMT:
      return collectMetaKeyFromShowLicence(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_VGROUPS_STMT:
      return collectMetaKeyFromShowVgroups(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_TOPICS_STMT:
      return collectMetaKeyFromShowTopics(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_SHOW_TRANSACTIONS_STMT:
      return collectMetaKeyFromShowTransactions(pCxt, (SShowStmt*)pStmt);
    case QUERY_NODE_DELETE_STMT:
      return collectMetaKeyFromDelete(pCxt, (SDeleteStmt*)pStmt);
    default:
      break;
  }
  return TSDB_CODE_SUCCESS;
}

int32_t collectMetaKey(SParseContext* pParseCxt, SQuery* pQuery, SParseMetaCache* pMetaCache) {
  SCollectMetaKeyCxt cxt = {.pParseCxt = pParseCxt, .pMetaCache = pMetaCache, .pStmt = pQuery->pRoot};
  return collectMetaKeyFromQuery(&cxt, pQuery->pRoot);
}