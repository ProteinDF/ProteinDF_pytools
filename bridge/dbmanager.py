#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import types
import logging
import sqlite3

class DbManager(object):
    """
    >>> db = DbManager()
    >>> db.create_table('table1', ['id', 'coulumn1', 'coulumn2'], 'id')
    >>> db.get_table_names()
    [u'table1']
    >>> db.has_table('table1')
    True

    >>> db.get_field_names('table1')
    ['id', 'coulumn1', 'coulumn2']
    >>> db.insert('table1', {'id':1, 'coulumn1':'Aichi', 'coulumn2':'Nagoya'})
    >>> db.insert('table1', {'id':2, 'coulumn1':'Miyagi', 'coulumn2':'Sendai'})
    >>> db.insert('table1', {'id':3, 'coulumn1':'Tokyo', 'coulumn2':'Tokyo'})

    >>> db.select('table1', where='coulumn1 = "Tokyo"')
    [{'coulumn1': u'Tokyo', 'coulumn2': u'Tokyo', 'id': 3}]

    >>> db.update('table1', contents={'coulumn2':'Shinjuku'}, \
        where='coulumn1 = "Tokyo"')

    >>> db.select('table1', where='coulumn1 = "Tokyo"')
    [{'coulumn1': u'Tokyo', 'coulumn2': u'Shinjuku', 'id': 3}]

    >>> db.select('table1')
    [{'coulumn1': u'Aichi', 'coulumn2': u'Nagoya', 'id': 1},\
 {'coulumn1': u'Miyagi', 'coulumn2': u'Sendai', 'id': 2},\
 {'coulumn1': u'Tokyo', 'coulumn2': u'Shinjuku', 'id': 3}]

    >>> db.delete('table1', where='id = 1')
    >>> db.select('table1')
    [{'coulumn1': u'Miyagi', 'coulumn2': u'Sendai', 'id': 2},\
 {'coulumn1': u'Tokyo', 'coulumn2': u'Shinjuku', 'id': 3}]
    """
    def __init__(self, db=':memory:', sql_debugout=False):
        self._logger = logging.getLogger(__name__)
        self._connection = sqlite3.connect(db)
        self._cursor = self._connection.cursor()
        self._sql_debugout = sql_debugout

    def __del__(self):
        self._connection.close()

    def __getitem__(self, key):
        answer = None
        if (self.has_table(key)):
            answer = DbTable(db_manager=self,
                             table_name=key)
        return answer

    # table ==================================================================
    def create_table(self, table_name, field_names, primary_key=None):
        """
        TABLEを作成する

        table_name: テーブル名
        field_names: フィールド(カラム)名
        型を指定しない場合は、field_namesをリストにする。
        型を指定する場合は、field_namesは{名前:型}の辞書型にする。
        """
        if not self.has_table(table_name):
            fields = []
            if isinstance(field_names, types.DictType):
                for k, v in field_names:
                    fields.append('{name} {type}'.format(k, v))
                field_names = fields
            fields_str = ', '.join(field_names)

            primary_key_str = ''
            if primary_key:
                if isinstance(primary_key, types.StringType):
                    primary_key_str = ', PRIMARY KEY({0})'.format(primary_key)
                if isinstance(primary_key, types.ListType):
                    primary_key_str += ', PRIMARY KEY ('
                    primary_key_str += ', '.join(primary_key) + ')'

            sql = 'CREATE TABLE {table_str} ({fields_str} {primary_key_str});'
            sql = sql.format(table_str=table_name,
                             fields_str=fields_str,
                             primary_key_str=primary_key_str)
            self.execute(sql)
        else:
            sys.stderr.write('already exist table: %s\n' % (table_name))

    def get_table_names(self):
        """
        TABLEの名前をリストで返す
        """
        table_names = []
        sql = "SELECT name FROM sqlite_master WHERE type='table'"
        self.execute(sql)
        results = self._cursor.fetchall()
        for row in results:
            table_name = row[0]
            table_names.append(table_name)
        return table_names

    def has_table(self, table_name):
        """
        指定されたTABLEを保持しているかどうかを返す
        """
        table_names = self.get_table_names()
        return (table_name in table_names)

    # field ==================================================================
    def get_field_names(self, table_name, fields="*"):
        """
        指定されたTABLE内のフィールドをリストで返す
        """
        field_names = None
        if (self.has_table(table_name) == True):
            field_names = []
            sql = "SELECT {0} FROM {1} LIMIT 1;".format(fields, table_name)
            self.execute(sql)
            for col, field_description in enumerate(self._cursor.description):
                field_name = field_description[0]
                field_names.append(field_name)
        return field_names

    def insert(self, table, contents):
        """
        データレコードを追加

        contentsは、filedをキーとした辞書型
        """
        fields = []
        values = []
        for field, value in contents.items():
            fields.append(field)
            values.append(value)
        fields_str = ", ".join(fields)
        values_str = ", ".join('?' for v in values)
        sql = "INSERT INTO {table}({fields}) VALUES({values});"
        sql = sql.format(table=table,
                         fields=fields_str,
                         values=values_str)
        self.execute(sql, values)
        self._connection.commit()

    def update(self, table, contents, where):
        """
        データレコードを更新

        contents、whereは、filedをキーとした辞書型
        """
        parameters = []
        set_sections = []
        for field, value in contents.items():
            set_sections.append('%s=?' % (field))
            parameters.append(value)
        set_str = ', '.join(set_sections)

        where_str = ''
        if isinstance(where, types.StringType):
            where_str = 'WHERE ' + where
        elif isinstance(where, type.DictType):
            where_sections = []
            for key, value in where.items():
                where_sections.append('%s=?' % (key))
                parameters.append(value)
            where_str = 'WHERE ' + ', '.join(where_sections)
        sql = 'UPDATE {table} SET {set_str} {where_str};'
        sql = sql.format(table=table,
                         set_str=set_str,
                         where_str=where_str)
        self.execute(sql, parameters)
        self._connection.commit()

    def delete(self, table, where):
        """
        データレコードを削除
        """
        parameters = []

        if isinstance(where, types.StringType):
            where_str = 'WHERE ' + where
        elif isinstance(where, types.DictType):
            where_sections = []
            for key, value in where.items():
                where_sections.append('%s=?' % (key))
                parameters.append(value)
            where_str = 'WHERE ' + ', '.join(where_sections)

        sql = 'DELETE FROM {table} {where_str};'
        sql = sql.format(table=table,
                         where_str=where_str)
        self.execute(sql, parameters)
        self._connection.commit()

    def select(self, table, fields=None, where=None):
        """
        データを取得する
        where句はANDのみサポート
        """
        # make SQL
        parameters = []

        field_str = '*'
        if fields:
            field_str = ', '.join(fields)

        where_str = ''
        if isinstance(where, types.StringType):
            where_str = 'WHERE ' + where
        elif isinstance(where, types.DictType):
            where_sections = []
            for key, value in where.items():
                where_sections.append('%s=?' % (key))
                parameters.append(value)
            where_str = 'WHERE ' + ' and '.join(where_sections)

        sql = 'SELECT {field_str} FROM {table} {where_str};'
        sql = sql.format(table=table,
                         field_str=field_str,
                         where_str=where_str)

        # execute
        self.execute(sql, parameters)

        data = self._cursor.fetchall()
        answer = [{}] * len(data)
        fields = []
        if self._cursor.description:
            for col, field_description in enumerate(self._cursor.description):
                field_name = field_description[0]
                fields.append(field_name)
            for row_index, row_data in enumerate(data):
                entry = {}
                for col_index, item in enumerate(row_data):
                    field_name = fields[col_index]
                    entry[field_name] = item
                answer[row_index] = entry
        return answer

    # SQL ====================================================================
    def execute(self, sql, parameters=None):
        """
        SQLを実行する
        """
        self._logger.debug("sql> {0}\n".format(sql))
        if (parameters != None):
            return self._cursor.execute(sql, parameters)
        else:
            return self._cursor.execute(sql)

    def get_results(self, sql):
        """
        SQLを実行し、結果をリストで返す
        """
        self.execute(sql)
        data = self._cursor.fetchall()
        field_names = []
        answer = []
        if self._cursor.description:
            for col, field_description in enumerate(self._cursor.description):
                field_name = field_description[0]
                field_names.append(field_name)

            for row in data:
                row_items = {}
                for index, item in enumerate(row):
                    row_items[field_names[index]] = item

            answer.append(row_items)
        return answer

    # etc ====================================================================
    def set_user_version(self, version):
        """
        ユーザーバージョンを設定する
        """
        version = int(version)
        self.execute('PRAGMA user_version = %d;' % (version))

    def get_user_version(self):
        """
        ユーザーバージョンを返す
        """
        answer = self.get_results('PRAGMA user_version;')
        answer = int(answer[0]['user_version'])
        return answer
        
    # output =================================================================
    def __str__(self):
        answer = ''
        tables = self.get_table_names()
        for table in tables:
            answer += self.pp_table(table)
        return answer

    def pp_table(self, table_name):
        """
        pretty print for table
        """
        answer = ''
        sql = "SELECT * FROM {0};".format(table_name)
        self.execute(sql)
        answer += 'TABLE: %s\n' % (table_name)
        answer += self.pp()
        #answer += '\n'
        return answer

    def pp(self, data=None, check_row_lengths=True):
        """
        pretty print for cursor data
        """
        if not data:
            data = self._cursor.fetchall()
        names = []
        lengths = []
        rules = []
        answer = ""
        if self._cursor.description:
            for col, field_description in enumerate(self._cursor.description):
                #print(field_description)
                field_name = field_description[0]
                names.append(field_name)
                field_length = field_description[2] or 12
                field_length = max(field_length, len(field_name))
                if check_row_lengths:
                    data_length = max([len(str(row[col])) for row in data])
                    field_length = max(field_length, data_length)
                lengths.append(field_length)
                rules.append('-' * field_length)
            format = " ".join(["%%-%ss" % l for l in lengths])
            result = [format % tuple(names), format % tuple(rules)]
            for row in data:
                result.append(format % tuple(row))
            answer = "\n".join(result)
        return answer

if __name__ == '__main__':
    import doctest
    doctest.testmod()
